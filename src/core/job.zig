const std = @import("std");
const math = std.math;

const print = std.debug.print;
const trace = @import("../tracy.zig").trace;

pub const Job = struct {
    pub const Error = error{JobError};
    pub const Result = enum { Complete, Retry, Abort };

    pub const ExecuteFunc = *const fn (self: *Job) Error!Result;
    pub const AbortFunc = *const fn (self: *Job) Error!void;

    executeFn: ExecuteFunc,
    abortFn: AbortFunc,

    // priority: i16,
    // name: []const u8,

    next: ?*Job = null,

    pub fn execute(self: *Job) Error!Result {
        // const ta = trace(@src());
        // defer ta.end();
        return try self.executeFn(self);
    }

    pub fn abort(self: *Job) Error!void {
        try self.abortFn(self);
    }

    pub fn implementor(self: *Job, comptime T: type, comptime field: []const u8) *T {
        return @fieldParentPtr(T, field, self);
    }

    pub fn implement(comptime T: type, e: ExecuteFunc, a: AbortFunc) T {
        return T{
            .job = Job{
                .executeFn = e,
                .abortFn = a,
            },
        };
    }
};

pub fn InPlaceQueue(comptime T: type) type {
    return struct {
        const Self = @This();

        head: ?*T = null,
        tail: ?*T = null,
        // lock:std.Thread.RwLock,

        //lock:std.Mutex,
        writeLock: std.Thread.Mutex,
        readLock: std.Thread.Mutex,
        wait: std.Thread.Semaphore,

        size: usize,

        pub fn init() Self {
            var s = Self{
                .head = null,
                .tail = null,
                .writeLock = std.Thread.Mutex{},
                .readLock = std.Thread.Mutex{},
                .wait = undefined,
                .size = 0,
            };

            return s;
        }

        pub fn tryPush(self: *Self, node: ?*T) bool {
            const ta = trace(@src());
            defer ta.end();

            if (!self.tryLockForWrite())
                return false;

            self.pushNoLock(node);
            self.unlockForWrite();
            return true;
        }

        pub fn tryLockForWrite(self: *Self) bool {
            if (!self.writeLock.tryLock())
                return false;

            self.readLock.lock();
            return true;
        }

        pub fn lockForWrite(self: *Self) void {
            self.writeLock.lock();

            self.readLock.lock();
        }

        pub fn unlockForWrite(self: *Self) void {
            self.writeLock.unlock();
            self.readLock.unlock();
        }
        pub fn pushNoLock(self: *Self, node: ?*T) void {
            const ta = trace(@src());
            defer ta.end();

            // print("push start {d}, {d}, {d}\n", .{@ptrToInt(self.head), @ptrToInt(self.tail), self.size});

            if (self.tail != null)
                self.tail.?.next = node;

            self.tail = node;

            _ = @atomicRmw(@TypeOf(self.size), &self.size, .Add, 1, .SeqCst);

            if (self.head == null) {
                self.head = self.tail;
                self.wait.post();
            }

            // if(self.wait.permits > 0)
            //     // print("push {d}, {d}, {d}\n", .{@ptrToInt(self.head), @ptrToInt(self.tail), self.size});
            //     self.wait.post();
        }

        pub fn push(self: *Self, node: ?*T) void {
            const ta = trace(@src());
            defer ta.end();
            self.lockForWrite();
            self.pushNoLock(node);
            self.unlockForWrite();
        }

        pub fn pop(self: *Self) ?*T {
            self.readLock.lock();
            defer self.readLock.unlock();

            var node = self.head orelse return null;
            if (self.head == self.tail) {
                self.head = null;
                self.tail = null;
            } else {
                self.head = self.head.?.next;
            }

            node.next = null;

            _ = @atomicRmw(@TypeOf(self.size), &self.size, .Sub, 1, .SeqCst);

            // if(self.head == null)
            //   self.wait.reset();

            return node;
        }

        pub fn popWait(self: *Self) ?*T {

            // const t = trace(@src());
            // defer t.end();

            var attempts: u16 = 0;
            while (true) {
                // const ta = trace(@src());
                // defer ta.end();

                var job = self.pop();
                attempts += 1;

                if (job != null) {
                    job.?.next = null;
                    return job;
                }
                // only wait if the queue is empty
                if (self.head == null and attempts > 63) {
                    // const tf = trace(@src());
                    // defer tf.end();
                    // print("pop wait {}, {d}\n", .{self.head, self.size});
                    self.wait.wait();
                    attempts = 0;
                }

                std.atomic.spinLoopHint();
            }
        }

        pub fn isEmpty(self: *Self) bool {
            return self.head == null;
        }

        pub fn count(self: *Self) usize {
            return @atomicLoad(@TypeOf(self.size), &self.size, .SeqCst);
        }
    };
}

pub const Queue = InPlaceQueue(Job);

//
pub const Worker = struct {
    pool: *WorkerPool,
    thread: std.Thread,
    id: u8,
    active: bool,
    pending: Queue,
    complete: Queue,
    executing: bool,

    pub fn init(ident: u8, pool: *WorkerPool) !Worker {
        var worker = Worker{
            .id = ident,
            .thread = undefined,
            .active = false,
            .pool = pool,
            .pending = Queue.init(),
            .complete = Queue.init(),
            .executing = false,
        };

        return worker;
    }

    pub fn start(self: *Worker) !void {
        self.thread = try std.Thread.spawn(.{}, Worker.run, .{self});
    }

    pub fn run(self: *Worker) void {
        self.active = true;
        while (self.active) {
            var job = self.pending.popWait();

            //self.runner.addRunningJob(job);

            @atomicStore(@TypeOf(self.executing), &self.executing, true, .SeqCst);

            const result = job.?.execute() catch continue;

            //self.complete.pushNoLock(job);

            @atomicStore(@TypeOf(self.executing), &self.executing, false, .SeqCst);

            //self.runner.removeRunningJob(job);
            _ = result;

            // switch (result) {
            //      .Complete, .Abort => self.complete.pushNoLock(job),
            //      //.Retry => self.runner.pending.push(job.?),
            //  }
        }
    }

    pub fn stop(self: *Worker) !void {
        self.active = false;
    }
};

// Group of workers pulling from the same queue
pub const WorkerPool = struct {
    workers: std.ArrayList(Worker),
    nextId: usize,

    pub const Counts = struct {
        pub fn init() Counts {
            return .{ .pending = 0, .running = 0, .complete = 0 };
        }

        pending: usize,
        running: usize,
        complete: usize,
    };

    pub fn init(allocator: *std.mem.Allocator, workerCount: u8) !WorkerPool {
        var workers = try std.ArrayList(Worker).initCapacity(allocator.*, workerCount);

        var self = WorkerPool{
            .workers = workers,
            .nextId = 0,
        };

        var w = workerCount;
        while (w > 0) {
            defer w -= 1;
            try self.workers.append(try Worker.init(w, &self));
        }

        return self;
    }

    /// Begin all jobs executing off of the queue
    pub fn start(self: *WorkerPool) !void {
        for (self.workers.items) |*w| {
            try w.start();
        }
    }

    pub fn stop(self: *WorkerPool) !void {
        for (self.workers.items) |*w| {
            try w.stop();
        }
    }

    pub fn deinit(self: *WorkerPool) void {
        self.stop();
        self.workers.deinit();
    }

    pub fn counts(self: *WorkerPool) Counts {
        var ret = Counts.init();
        for (self.workers.items) |*w| {
            //ret.complete += w.complete.count();
            ret.running += @intFromBool(w.executing);
            ret.pending += w.pending.count();
        }

        return ret;
    }

    pub fn push(self: *WorkerPool, job: *Job) void {
        // round robin queue
        self.nextId = (self.nextId + 1) % self.workers.items.len;

        while (!self.workers.items[self.nextId].pending.tryPush(job)) {
            self.nextId = (self.nextId + 1) % self.workers.items.len;
        }
    }
};

// pub fn Batch(comptime TDataType: type) type {
//     return struct {
//         const Self = @This();

//         pub const Error = error{BatchError};
//         pub const Result = enum { Complete, Retry, Abort };
//         pub const ExecuteFunc =*const fn (index:usize, total:usize, self: *TDataType) Error!Result;

//         work: std.ArrayList(TDataType),

//         executeFn: ExecuteFunc,

//         iterator: usize,
//         complete: usize,

//         pub fn init(executeFn: ExecuteFunc, allocator: std.mem.Allocator, capacity: usize) !Self {
//             var list = try std.ArrayList(TJobType).initCapacity(allocator, size);

//             var self = Self{
//                 .work = list,
//                 .iterator = 0,
//                 .complete = 0,
//                 .executeFn = executeFn,
//             };

//             return self;
//         }

//         pub fn reset(self: *Self) void {
//             // TODO: make and return error?
//             assert(!self.hasWork() and self.isComplete());

//             _ = @atomicStore(@TypeOf(self.iterator), &self.iterator, 0, .SeqCst);
//             _ = @atomicStore(@TypeOf(self.complete), &self.complete, 0, .SeqCst);
//         }

//         pub fn get(self: *Self, index: usize) ?*TJobType {
//             return &self.work.items[index] orelse null;
//         }

//         pub fn size(self: Self) usize {
//             return self.work.items.len;
//         }

//         pub fn resize(self: Self) !void {
//             try self.work.resize(size);
//         }

//         pub fn workerRunOne(self: *Self) !void {
//             const index = @atomicRmw(usize, &self.iterator, .Add, 1, .SeqCst);
//             var work = self.get(index);
//             if (work == null)
//                 return;

//             _ = try self.executeFn(index, self.work.items.len, work);

//             _ = @atomicRmw(usize, &self.complete, .Add, 1, .SeqCst);
//         }

//         pub fn isComplete(self: Self) bool {
//             return @atomicLoad(@TypeOf(self.complete), &self.complete, .SeqCst) >= self.size();
//         }

//         pub fn hasWork(self:Self) bool {
//             return @atomicLoad(@TypeOf(self.complete), &self.complete, .SeqCst) < self.size();
//         }
//     };
// }

pub const Runner = struct {
    const Self = @This();

    pending: Queue,
    running: usize,

    pub fn init() Self {
        return Self{ .pending = Queue.init(), .running = 0 };
    }

    pub fn addRunningJob(self: *Self, job: ?*Job) void {
        _ = job;
        _ = @atomicRmw(@TypeOf(self.running), &self.running, .Add, 1, .SeqCst);
    }

    pub fn removeRunningJob(self: *Self, job: ?*Job) void {
        _ = job;
        _ = @atomicRmw(@TypeOf(self.running), &self.running, .Sub, 1, .SeqCst);
    }

    pub fn count(self: *Self) usize {
        return @atomicLoad(@TypeOf(self.running), &self.running, .SeqCst) + self.pending.count();
    }
};

pub fn Pool(comptime TItemType: type) type {
    return struct {
        const Self = @This();
        const ItemList = std.ArrayList(TItemType);
        const AtomicUSize = std.atomic.Atomic(usize);

        items: ItemList,

        // nextItem:std.atomic.Atomic(usize),
        nextItem: usize,

        pub fn init(alloc: std.mem.Allocator, size: usize) !Self {
            var self = Self{
                // .nextItem = AtomicUSize.init(0),
                .nextItem = 0,
                .items = try ItemList.initCapacity(alloc, size),
            };

            try self.items.resize(size);

            return self;
        }

        pub fn getItem(self: *Self) *TItemType {
            const id = @atomicRmw(@TypeOf(self.nextItem), &self.nextItem, .Add, 1, .SeqCst);
            return &self.items.items[id];
        }

        pub fn reset(self: *Self) void {
            @atomicStore(@TypeOf(self.nextItem), &self.nextItem, 0, .SeqCst);
            // self.nextItem.store(0, .SeqCst);
        }

        pub fn deinit(self: *Self) void {
            self.items.deinit();
        }
    };
}

// ============================================================================
// Tests
// ============================================================================
pub const TestJob = struct {
    job: Job,
    id: u8,

    pub fn init(ident: u8) TestJob {
        return TestJob{
            .job = Job{
                // .priority = 0,
                // .name = "testjob",
                .executeFn = execute,
                .abortFn = abort,
            },

            .id = ident,
        };
    }

    fn execute(job: *Job) !Job.Result {
        const self = @fieldParentPtr(TestJob, "job", job);
        std.debug.print("\t job: {} execution!\n", .{self.id});
        return Job.Result.Complete;
    }

    fn abort(job: *Job) !void {
        const self = @fieldParentPtr(TestJob, "job", job);
        _ = self;
    }
};

pub const TestJob2 = struct {
    job: Job,

    pub fn init() TestJob2 {
        return Job.implement(TestJob2, execute, abort);
    }

    fn execute(job: *Job) !Job.Result {
        _ = job;
        // const self = @fieldParentPtr(TestJob2, "job", job);
        print("execution2!\n", .{});
        return Job.Result.Complete;
    }

    fn abort(job: *Job) !void {
        _ = job;
        // const self = @fieldParentPtr(TestJob2, "job", job);
    }
};

const assert = @import("std").debug.assert;
test "Job Interface Basic Execute" {
    var testjob = TestJob.init(0);
    var job: *Job = &testjob.job;
    _ = try job.execute();
    try job.abort();

    var testjob2 = TestJob2.init();
    job = &testjob2.job;
    _ = try job.execute();
    try job.abort();
}

fn anonExec(self: *Job) !Job.Result {
    std.debug.print("execution anon exec!\n", .{});
    _ = self;
    return Job.Result.Complete;
}

fn anonAbort(self: *Job) !void {
    _ = self;
    std.debug.print("execution anon abort!\n", .{});
}

test "Anonymous Job Interface Basic" {
    var j = Job{
        // .priority = 0,
        // .name = "Anon Job 1",
        .executeFn = anonExec,
        .abortFn = anonAbort,
    };

    var job: *Job = &j;

    _ = try job.execute();
    try job.abort();

    var testjob2 = TestJob2.init();
    job = &testjob2.job;
    _ = try job.execute();
    try job.abort();
}

test "Job Queue" {
    var j = Job{
        // .priority = 0,
        // .name = "Anon Job 1",
        .executeFn = anonExec,
        .abortFn = anonAbort,
    };

    var queue = InPlaceQueue(Job).init();

    var job: *Job = &j;

    queue.push(job);
    var pj = queue.pop();

    assert(pj == job);
    _ = try pj.?.execute();
}

test "Job Workers" {
    var queue = Runner.init();

    var workers = [_]Worker{
        try Worker.init(0, &queue),
        try Worker.init(1, &queue),
        try Worker.init(2, &queue),
        try Worker.init(3, &queue),
    };

    var jobs = [_]TestJob{
        TestJob.init(0),
        TestJob.init(1),
        TestJob.init(2),
        TestJob.init(3),
        TestJob.init(4),
        TestJob.init(5),
        TestJob.init(6),
        TestJob.init(7),
    };

    std.time.sleep(1_000_000);

    for (jobs) |*j| {
        queue.pending.push(&j.job);
    }

    for (workers) |*w| {
        try w.start();
    }

    while (queue.count() > 0) {
        // std.debug.print("remaining: {d}\n", .{queue.count()});
        std.time.sleep(1_000_000);
    }

    //std.os.sleep(100);
    assert(queue.count() == 0);
}

test "Job Worker Pool" {
    var allocator = std.heap.page_allocator;
    var queue = Runner.init();
    var pool = try WorkerPool.init(&queue, &allocator, 8);

    var jobs = [_]TestJob{
        TestJob.init(0),
        TestJob.init(1),
        TestJob.init(2),
        TestJob.init(3),
        TestJob.init(4),
        TestJob.init(5),
        TestJob.init(6),
        TestJob.init(7),
    };

    std.time.sleep(1_000_000);

    for (jobs) |*j| {
        queue.pending.push(&j.job);
    }

    try pool.start();

    while (queue.count() > 0) {
        // print("remaining: {d}\n", .{queue.count()});
        std.time.sleep(1_000_000);
    }

    try pool.stop();
    //std.os.sleep(100);
    assert(queue.count() == 0);
}
