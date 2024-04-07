const std = @import("std");
const math = std.math;
const assert = std.debug.assert;
const Thread = std.Thread;
const Timer = std.time.Timer;
const Instant = std.time.Instant;
const File = std.fs.File;
const json = std.json;

// TODO: make these parameters compile time, like generics
const SamplePoolCount = 2048;
const ThreadPoolCount = 32;

/// Track wall clock timing for sub steps (functions) of a single thread
pub const Profile = struct {
    nextSample: u32,
    depth: u8,
    samples: std.ArrayList(Sample),
    frameCount: u32,
    frameStartTime: Instant,

    /// Single timing block
    pub const Sample = struct {
        depth: u8,
        tag: []const u8,
        begin: Instant,
        end: Instant,
    };

    pub fn init(allocator: *std.mem.Allocator) !Profile {
        var list = try std.ArrayList(Sample).initCapacity(allocator.*, SamplePoolCount);
        return Profile{
            .nextSample = 1,
            .depth = 0,
            .frameCount = 0,
            .frameStartTime = undefined,
            .samples = list, // [_]Profile.Sample{ .{.depth=0, .tag="", .begin=0, .end=0} }** SamplePoolCount,
        };
    }

    /// Start tracking wall clock time
    pub fn beginSample(self: *Profile, tag: []const u8) u32 {
        const id = @atomicRmw(u32, &self.nextSample, .Add, 1, .SeqCst);

        var sample = self.samples.addOne() catch return 0;
        assert(id < self.samples.items.len);

        sample.depth = @atomicRmw(u8, &self.depth, .Add, 1, .SeqCst);
        sample.tag = tag;
        //TODO: find out how to correctly handle a timer error here

        sample.begin = Instant.now() catch return id;
        return id;
    }

    pub fn sampleTime(self: *Profile, id: u32) i64 {
        return @as(i64, @intCast(self.samples.items[id].end - self.samples.items[id].begin));
    }

    // Stop tracking wall clock timing
    pub fn endSample(self: *Profile, id: u32) void {
        // const d = @atomicRmw(u8, &self.depth, .Sub, 1, .SeqCst);

        //TODO: find out how to correctly handle a timer error here
        self.samples.items[id].end = Instant.now() catch return;
    }

    /// Reset profile data for a new frame
    pub fn nextFrame(self: *Profile) void {
        self.depth = 0;
        self.nextSample = 1;
        self.frameCount += 1;
        self.samples.resize(1) catch return;

        //TODO: find out how to correctly handle a timer error here
        self.frameStartTime = Instant.now() catch return;
    }

    pub fn hasSamples(self: Profile) bool {
        return self.nextSample > 1;
    }

    pub fn streamPrint(self: *Profile, file: File) !void {
        var stream = file.outStream();
        try stream.print("f:{}, sc:{}\n", .{ self.frameCount, self.nextSample });

        for (self.samples.items, 0..) |sample, i| {
            if (i > self.nextSample)
                break;

            if (i == 0 or sample.begin.timestamp == 0)
                continue;

            var s = sample.depth;
            while (s > 0) {
                try stream.print(" ", .{});
                s -= 1;
            }

            const begin = sample.begin.since(self.frameStartTime);
            const end = sample.end.since(self.frameStartTime);
            try stream.print("[{}:{}] b:{} ns, e:{} ns, d:{} ns, t:{}\n", .{ i, sample.depth, begin, end, end - begin, sample.tag });
        }
    }

    pub fn jsonPrint(self: *Profile, file: File) !void {
        var stream = file.writer();
        try stream.print("{{ \"traceEvents\":[", .{});

        for (self.samples.items, 0..) |sample, i| {
            if (i > self.nextSample)
                break;

            const order = sample.begin.order(self.frameStartTime);
            if (i == 0 or order == std.math.Order.lt)
                continue;

            const begin = sample.begin.since(self.frameStartTime);
            const end = sample.end.since(self.frameStartTime);

            try stream.print("{{ \"name\": \"{s}\", \"cat\":\"{s}\",\"ph\":\"{s}\",\"pid\": {any},\"tid\": {any},\"ts\": {any}, \"dur\":{any} }}", .{ sample.tag, "PERF", "X", 1, 1, begin, end - begin });

            //std.debug.warn("{}, {}\n", .{i, self.nextSample});
            if (i < self.nextSample - 1)
                try stream.print(",\n", .{});
        }

        try stream.print("]}}", .{});
    }

    pub fn jsonFileWrite(self: *Profile, allocator: std.mem.Allocator, filepath: []const u8) !void {
        const cwd = std.fs.cwd();

        var resolvedPath = try std.fs.path.resolve(allocator, &[_][]const u8{filepath});
        defer allocator.free(resolvedPath);

        std.debug.print("path: {s}", .{resolvedPath});

        var file = try cwd.createFile(resolvedPath, .{});
        defer file.close();

        try self.jsonPrint(file);
    }

    pub fn dumpStreamText(self: *Profile, file: File) !void {
        var stream = file.outStream();
        try stream.print("Frame, Frame Start, Id, Depth, Begin ns, End ns, Exec ns, tag\n", .{});

        for (self.samples.items, 0..) |sample, i| {
            if (i == 0 or sample.begin == 0)
                continue;

            const begin = sample.begin - self.frameStartTime;
            const end = sample.end - self.frameStartTime;
            try stream.print("{any}, {any}, {any}, {any}, {any}, {any}, {any}, {any}\n", .{ self.frameCount, self.frameStartTime, i, sample.depth, sample.begin, sample.end, end - begin, sample.tag });
        }
    }
};

pub const Sampler = struct {
    id: u32,
    owner: *Profile,
    tag: []const u8,

    pub fn init(profiler: *Profile, tag: []const u8) Sampler {
        return Sampler{ .owner = profiler, .id = 0, .tag = tag };
    }

    pub fn initAndBegin(profiler: *Profile, tag: []const u8) Sampler {
        return Sampler{
            .owner = profiler,
            .tag = tag,
            .id = profiler.beginSample(tag),
        };
    }

    pub fn begin(self: *Sampler) void {
        self.id = self.profiler.beginSample(self.tag);
    }

    pub fn end(self: *Sampler) void {
        self.owner.endSample(self.id);
    }

    pub fn value(self: *Sampler) i64 {
        return self.owner.sampleTime(self.id);
    }
};

test "Profiler.init" {
    // var profiler = ThreadProfile{
    //       .id = Thread.getCurrentId(),
    //       .nextSample = 1,
    //       .depth = 0,
    //       .samples = [_]ThreadProfile.Sample{ .{.depth=0, .tag="", .begin=0, .end=0} }** SamplePoolCount,
    //     };
    // {
    //   var sampleA = profiler.beginSample("A");
    //   defer profiler.endSample(sampleA);

    //   var i:u16 = 0;
    //   while(i < 10) {
    //     var sampleB = profiler.beginSample("B");
    //     defer profiler.endSample(sampleB);
    //     i += 1;
    //   }
    // }

    // std.debug.warn("\n",.{});
    // profiler.print();
}
