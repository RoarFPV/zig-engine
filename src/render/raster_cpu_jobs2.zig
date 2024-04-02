// compile in ubuntu:
// $ zig build-exe paint.zig --library SDL2 --library SDL2main --library c -isystem "/usr/include" --library-path "/usr/lib/x86_64-linux-gnu"

const std = @import("std");
const warn = std.debug.print;
const fmt = std.fmt;
const assert = @import("std").debug.assert;
const math = std.math;

const Vec4f = @import("../core/vector.zig").Vec4f;
const Mat44f = @import("../core/matrix.zig").Mat44f;

const Profile = @import("../core/profiler.zig").Profile;
pub const Font = @import("font.zig").Font;

const Mesh = @import("mesh.zig").Mesh;
const pixelbuffer = @import("pixel_buffer.zig");
const PixelBuffer = pixelbuffer.PixelBuffer;
const PixelRenderer = pixelbuffer.PixelRenderer;

const jobs = @import("../core/job.zig");
const Job = jobs.Job;
const JobWorker = jobs.Worker;
const JobQueue = jobs.Queue;
const JobRunner = jobs.Runner;
const JobPool = jobs.Pool;

const sdf = @import("../core/sdf.zig");

pub const material = @import("material.zig");
pub const Material = material.Material;

const tracy = @import("../tracy.zig");
const trace = tracy.trace;

/// RGBA 32 bit color value
pub const Color = struct {
    color: [4]u8 = [4]u8{ 0, 0, 0, 0 },

    pub fn r(self: Color) u8 {
        return self.color[0];
    }
    pub fn g(self: Color) u8 {
        return self.color[1];
    }
    pub fn b(self: Color) u8 {
        return self.color[2];
    }
    pub fn a(self: Color) u8 {
        return self.color[3];
    }
    pub fn setR(self: *Color, val: u8) void {
        self.color[0] = val;
    }
    pub fn setG(self: *Color, val: u8) void {
        self.color[1] = val;
    }
    pub fn setB(self: *Color, val: u8) void {
        self.color[2] = val;
    }
    pub fn setA(self: *Color, val: u8) void {
        self.color[3] = val;
    }

    pub fn white() Color {
        var color = Color{ .color = [4]u8{ 255, 255, 255, 255 } };
        return color;
    }

    pub fn black() Color {
        var color = Color{ .color = [4]u8{ 0, 0, 0, 255 } };
        return color;
    }

    pub fn init(cr: u8, cg: u8, cb: u8, ca: u8) Color {
        return Color{ .color = [4]u8{ cr, cg, cb, ca } };
    }

    inline fn norm(e: f32) u8 {
        return @as(u8, @intFromFloat(@min(1.0, @max(0.0, e)) * 255));
    }

    pub inline fn fromNormal(cr: f32, cg: f32, cb: f32, ca: f32) Color {
        return Color.init(norm(cr), norm(cg), norm(cb), norm(ca));
    }

    pub inline fn fromNormalVec4f(vec: Vec4f) Color {
        return Color.fromNormal(vec.x(), vec.y(), vec.z(), vec.w());
    }
};

const Bounds = struct {
    min: Vec4f,
    max: Vec4f,

    pub fn init(min: Vec4f, max: Vec4f) Bounds {
        return Bounds{
            .min = min,
            .max = max,
        };
    }

    pub fn add(self: *Bounds, point: Vec4f) void {
        self.min = Vec4f.min(self.min, point);
        self.max = Vec4f.max(self.max, point);
    }

    pub fn limit(self: *Bounds, l: Bounds) void {
        self.min = Vec4f.max(l, self.min);
        self.max = Vec4f.min(l, self.max);
    }

    pub fn topLeftHandLimit(self: *Bounds) void {
        self.min.sub(Vec4f.half());
        self.min.ceil();

        self.max.sub(Vec4f.half());
        self.max.ceil();
    }
};

const ColorPixelBuffer = PixelBuffer(Color);
const DepthPixelBuffer = PixelBuffer(f32);

var colorBuffer: [2]ColorPixelBuffer = undefined;
var depthBuffer: [2]DepthPixelBuffer = undefined;

var currentBuffer: u1 = 0;

var profile: ?*Profile = undefined;
var allocator: *std.mem.Allocator = undefined;
var colorRenderer: PixelRenderer(Color) = undefined;

pub const Stats = struct {
    const Self = @This();

    queuedMeshes: u32 = 0,
    totalMeshes: u32 = 0,
    totalTris: u32 = 0,
    totalPixels: u32 = 0,

    renderedMeshes: u32 = 0,
    renderedTris: u32 = 0,
    renderedPixels: u32 = 0,

    trisTooSmall: u32 = 0,
    trisTooNear: u32 = 0,
    trisBackfacing: u32 = 0,
    jobWaitCount: u32 = 0,
    jobCount: u32 = 0,

    pub fn init() Self {
        return Self{
            .totalMeshes = 0,
            .totalTris = 0,
            .totalPixels = 0,
            .renderedMeshes = 0,
            .renderedTris = 0,
            .renderedPixels = 0,
            .trisTooSmall = 0,
            .trisTooNear = 0,
            .trisBackfacing = 0,
            .jobWaitCount = 0,
            .jobCount = 0,
            .queuedMeshes = 0,
        };
    }

    pub fn reset(self: *Self) void {
        const ti = @typeInfo(Self);
        inline for (ti.Struct.fields) |field| {
            @field(self, field.name) = 0;
        }
    }

    pub fn trace(self: Stats) void {
        const ti = @typeInfo(Stats);
        inline for (ti.Struct.fields) |field| {
            tracy.plotValue(field.name.ptr, @field(self, field.name));
        }
    }

    pub fn print(self: *Self) void {
        std.debug.print(" [ m({}, {}|{d:.2}%), t({}, {}|{d:.2}%, <{d}, |<{d}, bf{d}), p({}, {}|{d:.2}%) ]\n", .{
            self.renderedMeshes,
            self.totalMeshes,
            @as(f32, @floatFromInt(self.renderedMeshes)) / @as(f32, @floatFromInt(self.totalMeshes)) * 100.0,

            self.renderedTris,
            self.totalTris,
            @as(f32, @floatFromInt(self.renderedTris)) / @as(f32, @floatFromInt(self.totalTris)) * 100.0,
            self.trisTooSmall,
            self.trisTooNear,
            self.trisBackfacing,

            self.renderedPixels,
            self.totalPixels,
            @as(f32, @floatFromInt(self.renderedPixels)) / @as(f32, @floatFromInt(self.totalPixels)) * 100.0,
        });
    }
};

pub const MeshRenderData = struct {
    model: *const Mat44f,
    view: *const Mat44f,
    proj: *const Mat44f,
    mv: Mat44f,
    vp: Mat44f,
    mvp: Mat44f,
    offset: u16,
    mesh: *Mesh,
    shader: *Material,
};

pub const TriRenderData = struct {
    pub const Visible = 1 << 0;
    pub const TooSmall = 1 << 1;

    id: u32,
    offset: u16,

    meshData: *MeshRenderData,

    indicies: [3]u16,

    color: [3]Vec4f,
    uv: [3]Vec4f,
    normals: [3]Vec4f,
    worldNormals: [3]Vec4f,

    rawVertex: [3]Vec4f,

    // transformed and projected verticies
    screenVertex: [3]Vec4f, // Projection * Model * View
    cameraVertex: [3]Vec4f, // Model * View
    worldVertex: [3]Vec4f, // Model

    normal: Vec4f,
    edges: [3]Vec4f,
    screenArea: f32,
    flags: u32,
    backfacing: bool,

    screenBounds: Bounds,
};

const TriSpanData = struct {
    triData: *TriRenderData,
};

// TODO: investigate using a variant for this
pub fn RenderJob(comptime TDataType: type) type {
    const RenderExecuteFunc = switch (@typeInfo(TDataType)) {
        .Pointer => *const fn (data: TDataType) void,
        else => *const fn (data: *TDataType) void,
    };

    return struct {
        const Self = @This();

        func: RenderExecuteFunc,
        job: Job,
        complete: bool = false,
        data: TDataType,

        pub fn init(func: RenderExecuteFunc) Self {
            return Self{
                .job = Job{
                    .executeFn = execute,
                    .abortFn = abort,
                    .next = null,
                },
                .func = func,
                .data = undefined,
                .complete = false,
            };
        }

        fn execute(job: *Job) Job.Error!Job.Result {
            const self = job.implementor(Self, "job");

            switch (@typeInfo(TDataType)) {
                .Pointer => self.func(self.data),
                else => self.func(&self.data),
            }

            self.complete = true;
            //std.debug.warn("\t job: {}:{} execution!\n", .{self.id, self.complete});
            return Job.Result.Complete;
        }

        fn abort(job: *Job) Job.Error!void {
            _ = job.implementor(Self, "job");
        }
    };
}

const MeshRenderJob = RenderJob(MeshRenderData);
const TriRenderJob = RenderJob(TriRenderData);
const SpanRenderJob = RenderJob(TriSpanData);

const ColorBufferClearJob = RenderJob(*ColorPixelBuffer);
const DepthBufferClearJob = RenderJob(*DepthPixelBuffer);

var meshJobs: JobPool(MeshRenderJob) = undefined;
var triJobs: JobPool(TriRenderJob) = undefined;
var spanJobs: JobPool(SpanRenderJob) = undefined;

// var renderQueue = JobRunner.init();
var renderWorkers: jobs.WorkerPool = undefined;
var stats = Stats.init();
var viewport = Vec4f.zero();
var renderBounds: Bounds = undefined;

var colorBufferClearJob: ColorBufferClearJob = undefined;
var depthBufferClearJob: DepthBufferClearJob = undefined;

pub fn frameStats() Stats {
    return stats;
}

pub fn drawLine(xFrom: i32, yFrom: i32, xTo: i32, yTo: i32, color: Color) void {
    const zone = trace(@src());
    defer zone.end();
    colorRenderer.drawLine(xFrom, yFrom, xTo, yTo, color);
}

pub fn bufferStart() *u8 {
    return &colorBuffer[currentBuffer].bufferStart().color[0];
}

pub fn bufferLineSize() usize {
    return colorBuffer[currentBuffer].bufferLineSize();
}

pub fn getViewport() Vec4f {
    return viewport;
}

pub fn init(renderWidth: u16, renderHeight: u16, alloc: *std.mem.Allocator, profileContext: ?*Profile) !void {
    profile = profileContext;
    allocator = alloc;

    const clearColor = Color.init(50, 50, 120, 255);
    const clearDepth = std.math.inf(f32);
    //colorBuffer.clear();
    //depthBuffer.clear(std.math.inf(f32));

    currentBuffer = 0;

    colorBuffer[0] = try PixelBuffer(Color).init(renderWidth, renderHeight, clearColor, allocator);
    depthBuffer[0] = try PixelBuffer(f32).init(renderWidth, renderHeight, clearDepth, allocator);

    colorBuffer[1] = try PixelBuffer(Color).init(renderWidth, renderHeight, clearColor, allocator);
    depthBuffer[1] = try PixelBuffer(f32).init(renderWidth, renderHeight, clearDepth, allocator);

    inline for (&colorBuffer) |*buffer| {
        buffer.clearDefault();
    }

    colorBufferClearJob = ColorBufferClearJob.init(ColorPixelBuffer.clearDefault);
    depthBufferClearJob = DepthBufferClearJob.init(DepthPixelBuffer.clearDefault);

    colorRenderer = PixelRenderer(Color).init(&colorBuffer[currentBuffer]);

    renderBounds = Bounds.init(Vec4f.zero(), viewport);
    viewport = Vec4f.init(@as(f32, @floatFromInt(colorBuffer[0].w)), @as(f32, @floatFromInt(colorBuffer[0].h)), 0, 0);

    meshJobs = try JobPool(MeshRenderJob).init(alloc.*, 16);
    triJobs = try JobPool(TriRenderJob).init(alloc.*, 1024 * 1024);
    spanJobs = try JobPool(SpanRenderJob).init(alloc.*, 1024);

    renderWorkers = try jobs.WorkerPool.init(alloc, @as(u8, @intCast(try std.Thread.getCpuCount())) / 2);
    try renderWorkers.start();
}

pub fn shutdown() void {
    colorBuffer[0].deinit();
    colorBuffer[1].deinit();
    depthBuffer[0].deinit();
    depthBuffer[1].deinit();
    meshJobs.deinit();
    triJobs.deinit();
    spanJobs.deinit();
}

pub fn beginFrame(profiler: ?*Profile) *u8 {
    const zone = trace(@src());
    defer zone.end();

    profile = profiler;
    var pprof = profile.?.beginSample("render.beginFrame");
    defer profile.?.endSample(pprof);

    stats.reset();

    colorBufferClearJob.data = &colorBuffer[currentBuffer];
    colorBufferClearJob.complete = false;

    depthBufferClearJob.data = &depthBuffer[currentBuffer];
    depthBufferClearJob.complete = false;

    renderWorkers.push(&colorBufferClearJob.job);
    renderWorkers.push(&depthBufferClearJob.job);

    currentBuffer = ~currentBuffer;

    colorRenderer.setBuffer(&colorBuffer[currentBuffer]);

    triJobs.reset();
    meshJobs.reset();
    spanJobs.reset();

    return &colorBuffer[currentBuffer].bufferStart().color[0];
}

pub fn endFrame() void {
    //pixels.swapBuffers();
    {
        const zone = trace(@src());
        defer zone.end();

        stats.jobWaitCount = 0;
        stats.jobWaitCount = 0;

        // if(renderQueue.count() > 0)
        // {
        //     while( stats.totalMeshes <= 0 )
        //     {
        //         //stats.jobCount+=1;
        //         std.atomic.spinLoopHint();
        //     }
        // }

        //const startCounts= renderWorkers.counts();

        while (stats.jobWaitCount < 1e6) {
            const counts = renderWorkers.counts();

            if (counts.pending <= 0 and counts.running <= 0)
                break;

            stats.jobWaitCount += 1;
            std.atomic.spinLoopHint();
        }

        //assert(stats.jobWaitCount > 0);
        //assert(stats.totalMeshes == stats.renderedMeshes);
    }

    //stats.trace();
}

pub fn drawMesh(m: *const Mat44f, v: *const Mat44f, p: *const Mat44f, mesh: *Mesh, shader: *Material) void {
    const zone = trace(@src());
    defer zone.end();

    var job = meshJobs.getItem();

    job.* = MeshRenderJob.init(drawMeshJob);

    job.data.mesh = mesh;
    job.data.model = m;
    job.data.view = v;
    job.data.proj = p;
    job.data.shader = shader;

    // Model * View
    var mv = v.*;
    mv.mul(m.*);
    job.data.mv = mv;

    // View * projection
    var vp = p.*;
    vp.mul(v.*);

    // Projection * Model * View
    var mvp = p.*;
    mvp.mul(v.*);
    mvp.mul(m.*);
    job.data.mvp = mvp;

    renderWorkers.push(&job.job);
}

fn drawMeshJob(meshJob: *MeshRenderData) void {
    const zone = trace(@src());
    defer zone.end();

    var t: u16 = 0;
    const ids = meshJob.mesh.indexBuffer.len;

    _ = @atomicRmw(u32, &stats.totalMeshes, .Add, 1, .SeqCst);

    // renderQueue.pending.lockForWrite();
    // defer renderQueue.pending.unlockForWrite();

    while (t < ids / 3) {
        var job = triJobs.getItem();
        job.* = TriRenderJob.init(drawTriJob);
        job.data.id = @as(u16, @truncate(t));
        job.data.offset = t * 3;
        job.data.meshData = meshJob;
        job.complete = false;
        renderWorkers.push(&job.job);
        t += 1;

        // if(t % 16 == 0)
        // {
        //     locked = false;
        //     renderQueue.pending.unlockForWrite();
        // }
    }

    _ = @atomicRmw(u32, &stats.renderedMeshes, .Add, 1, .SeqCst);
}

fn getTriangleNormal(points: [3]Vec4f) Vec4f {
    var e0 = points[1];
    var e1 = points[2];

    e0.sub(points[0]);
    e1.sub(points[0]);

    return e0.cross3(e1).normalized3();
}

fn drawTriJob(triJob: *TriRenderData) void {
    const zone = trace(@src());
    defer zone.end();

    var tri = triJob;
    const data = tri.meshData;
    const mesh = data.mesh;
    const shader = data.shader;
    // var sortedByY:[3]u8 = .{};

    _ = @atomicRmw(u32, &stats.totalTris, .Add, 1, .SeqCst);

    comptime var p = 0;
    inline while (p < 3) {
        var offset = tri.offset + p;
        assert(offset < mesh.indexBuffer.len);
        tri.indicies[p] = mesh.indexBuffer[offset];
        tri.rawVertex[p] = mesh.vertexBuffer[tri.indicies[p]];

        tri.cameraVertex[p] = data.mv.mul33_vec4(tri.rawVertex[p]);
        tri.screenVertex[p] = shader.projectionShader(data.proj, shader.vertexShader(&data.mv, offset, tri.rawVertex[p], shader), viewport, shader);

        tri.normals[p] = mesh.vertexNormalBuffer[mesh.indexNormalBuffer[offset]];
        tri.worldNormals[p] = data.model.mul33_vec4(tri.normals[p]);
        tri.uv[p] = Vec4f.init(mesh.textureCoordBuffer[mesh.indexUVBuffer[(offset)] * 2 + 0], mesh.textureCoordBuffer[mesh.indexUVBuffer[(offset)] * 2 + 1], 0, 0);

        tri.worldVertex[p] = data.model.mul_vec4(tri.rawVertex[p]);

        const indexScalar = @as(f32, @floatFromInt(offset)) / @as(f32, @floatFromInt(mesh.indexBuffer.len));
        tri.color[p] = Vec4f.init(0.4, 0.7, 0.5, 1.0).mulDup(Vec4f.init(indexScalar, indexScalar, indexScalar, 1.0)); //mesh.colorBuffer[tri.indicies[p]];

        p += 1;
    }

    tri.normal = getTriangleNormal(tri.cameraVertex);
    tri.backfacing = tri.cameraVertex[0].normalized3().dot3(tri.normal) > 0.0000000000001;
    tri.screenArea = Vec4f.triArea(tri.screenVertex[0], tri.screenVertex[1], tri.screenVertex[2]);
    tri.screenBounds = Bounds.init(tri.screenVertex[0], tri.screenVertex[0]);

    if (tri.backfacing) {
        _ = @atomicRmw(u32, &stats.trisBackfacing, .Add, 1, .SeqCst);
        return;
    }

    p = 1;
    inline while (p < 3) {
        tri.screenBounds.add(tri.screenVertex[p]);
        p += 1;
    }

    //tri.screenBounds.limit(renderBounds);
    tri.screenBounds.topLeftHandLimit();

    // Too small to see
    if (tri.screenArea <= 0) {
        _ = @atomicRmw(u32, &stats.trisTooSmall, .Add, 1, .SeqCst);
        return;
    }

    for (tri.screenVertex) |v| {
        if (v.z() <= 0.1) {
            _ = @atomicRmw(u32, &stats.trisTooNear, .Add, 1, .SeqCst);
            return;
        }
    }

    // drawLine(
    //     @floatToInt(i32, tri.screenVertex[0].x()),
    //     @floatToInt(i32, tri.screenVertex[0].y()),
    //     @floatToInt(i32, tri.screenVertex[1].x()),
    //     @floatToInt(i32, tri.screenVertex[1].y()), Color.fromNormalVec4f(tri.color[0]));

    // drawLine(
    //     @floatToInt(i32, tri.screenVertex[1].x()),
    //     @floatToInt(i32, tri.screenVertex[1].y()),
    //     @floatToInt(i32, tri.screenVertex[2].x()),
    //     @floatToInt(i32, tri.screenVertex[2].y()), Color.fromNormalVec4f(tri.color[1]));

    // drawLine(
    //     @floatToInt(i32, tri.screenVertex[2].x()),
    //     @floatToInt(i32, tri.screenVertex[2].y()),
    //     @floatToInt(i32, tri.screenVertex[0].x()),
    //     @floatToInt(i32, tri.screenVertex[0].y()), Color.fromNormalVec4f(tri.color[2]));

    // math.max(math.max(v[0].y(), v[1].y()), v[2].y());

    shadeTriJob(tri);

    _ = @atomicRmw(u32, &stats.renderedTris, .Add, 1, .SeqCst);
}

fn shadeTriJob(triData: *TriRenderData) void {
    const bounds = triData.screenBounds;
    var y = bounds.min.y();
    var p: Vec4f = Vec4f.init(0, 0, 0, 0);
    var pixelNormal: Vec4f = Vec4f.init(0, 0, 0, 0);
    var fbc: Vec4f = Vec4f.init(0, 0, 0, 1);
    var uv: Vec4f = Vec4f.init(0, 0, 0, 1);
    var c: Color = Color.black();

    _ = @atomicRmw(u32, &stats.renderedTris, .Add, 1, .SeqCst);

    // Generate spans
    while (y <= bounds.max.y()) {
        var x = bounds.min.x();
        defer y += 1;

        while (x <= bounds.max.x()) {
            // const pixtracy = trace(@src());
            // defer pixtracy.end();
            // var pprof = profile.?.beginSample("render.mesh.draw.tri.pixel");
            // defer profile.?.endSample(pprof);
            //_ = @atomicRmw(u32, &stats.totalPixels, .Add, 1, .SeqCst);

            defer x += 1;

            p.setX(x);
            p.setY(y);

            var triBary = Vec4f.triBarycentericCoordsOld(triData.screenVertex[0], triData.screenVertex[1], triData.screenVertex[2], p);

            // TODO: near plane clipping
            if (triBary.x() < 0 or triBary.y() < 0 or triBary.z() < 0)
                continue;

            triBary.div3(triData.screenArea);

            // if we use perspective correct interpolation we need to
            // multiply the result of this interpolation by z, the depth
            // of the point on the 3D triangle that the pixel overlaps.
            const z = (triBary.x() * triData.screenVertex[0].z() + triBary.y() * triData.screenVertex[1].z() + triBary.z() * triData.screenVertex[2].z());

            if (triData.meshData.shader.depthTest == 1 and
                depthBuffer[currentBuffer]
                .setLessThan(@as(i32, @intFromFloat(x)), @as(i32, @intFromFloat(y)), z) == 0)
                continue;

            p.setZ(z);

            // interpolate vertex colors across all pixels
            fbc = Vec4f.triInterpArray(triBary, triData.color, 1.0, 1.0);
            pixelNormal = Vec4f.triInterpArray(triBary, triData.worldNormals, 1.0, 1.0);
            uv = Vec4f.triInterpArray(triBary, triData.uv, 1.0, 1.0);

            var vc = triData.meshData.shader.pixelShader(&triData.meshData.mvp, p, fbc, pixelNormal, uv, triData.meshData.shader);

            if (vc.w() <= 0.0)
                continue;

            vc.clamp01();
            vc.scale(255);

            c.setR(@as(u8, @intFromFloat(@fabs(vc.x()))));
            c.setG(@as(u8, @intFromFloat(@fabs(vc.y()))));
            c.setB(@as(u8, @intFromFloat(@fabs(vc.z()))));
            c.setA(@as(u8, @intFromFloat(@fabs(vc.w()))));

            //_ = @atomicRmw(u32, &stats.renderedPixels, .Add, 1, .SeqCst);
            writePixel(@as(i32, @intFromFloat(x)), @as(i32, @intFromFloat(y)), z, c);
        }
    }
}

fn drawTriSpanJob(spanJob: *TriSpanData) void {
    _ = spanJob;
}

pub fn drawMesh_old(m: *const Mat44f, v: *const Mat44f, p: *const Mat44f, mesh: *Mesh, shader: *Material) void {
    // const tracy = trace(@src());
    // defer tracy.end();
    var sp = profile.?.beginSample("render.mesh.draw");
    defer profile.?.endSample(sp);

    _ = @atomicRmw(u32, &stats.totalMeshes, .Add, 1, .SeqCst);
    _ = shader;

    var mv = m.*;
    mv.mul(v.*);

    var mvp = mv.*;
    mvp.mul(p.*);
    // mvp.mul(m.*);

    const ids = mesh.indexBuffer.len;
    // const numTris = ids / 3;

    var t: u16 = 0;

    // while (t < ids) {
    //     triJobs.items[t / 3] = TriRenderJob.init(@truncate(u8, t / 3), m, v, p, &mv, &mvp, t, mesh, shader);

    //     var data = &triJobs.items[t / 3];
    //     data.complete = false;
    //     renderQueue.push(&data.job) catch continue;
    //     //drawTri(data);
    //     t += 3;
    // }

    t = 0;
    while (t < ids) {
        var wait: u64 = 0;
        var wt = profile.?.beginSample("render.mesh.wait.tri");
        defer profile.?.endSample(wt);
        var job: *TriRenderJob = &triJobs.items[t / 3];

        while (!job.complete) // and wait < 1_000_000)
        {
            std.atomic.spinLoopHint();
            // std.SpinLock.yield();
            //std.debug.warn("\t job: {}:{}:{} waiting!\n", .{job.id, job.complete, wait});
            wait += 1;
        }

        t += 3;
    }

    _ = @atomicRmw(u32, &stats.renderedMeshes, .Add, 1, .SeqCst);
}

pub fn drawPointMesh(mvp: *const Mat44f, mesh: *Mesh, shader: *Material) void {
    // const ids = mesh.vertexBuffer.len;
    _ = shader;
    for (mesh.vertexBuffer, 0..) |vertex, i| {
        drawPoint(mvp, mvp, vertex, mesh.colorBuffer[i]);
    }
}

pub fn drawString(font: *Font, str: []const u8, x: i32, y: i32, color: Vec4f) void {
    const zone = trace(@src());
    defer zone.end();
    const colorValue = Color.init(@as(u8, @intFromFloat(color.x() * 255)), @as(u8, @intFromFloat(color.y() * 255)), @as(u8, @intFromFloat(color.z() * 255)), @as(u8, @intFromFloat(color.w() * 255)));

    const message = str;
    _ = colorValue;
    // _ = shader;

    var lines: u8 = 0;
    var startx: i32 = x;
    var drawCount: u8 = 0;

    for (message) |c| {
        const cx = font.characterX(c);
        const cy = font.characterY(c);

        if (c == '\n') {
            lines += 1;
            startx = x;
            drawCount = 0;
            continue;
        }

        // TODO: copy lines rather than sample directly?

        var starty = y + lines * font.glyphHeight + 2;

        var coy: i32 = 1;

        var baseX = (x + (font.glyphWidth - 1) * @as(i32, @intCast(drawCount)));

        while (coy <= font.glyphHeight) {
            defer coy += 1;

            var cox: i32 = 1;
            while (cox <= font.glyphWidth + 1) {
                defer cox += 1;

                var samp = font.characterColor(cx, cy, cox, coy);
                if (samp.x() <= 0.001)
                    continue;

                writePixelNormal(baseX + cox, starty + coy, 1.0, samp);
            }
        }

        drawCount += 1;
    }
}

//
pub fn drawPoint(proj: *const Mat44f, mv: *const Mat44f, point: Vec4f, color: Vec4f, shader: *Material) void {
    const px = shader.projectionShader(proj, shader.vertexShader(mv, 0, point, shader), viewport, shader);

    //   shader.vertexShader(mvp, 0, point, shader);
    const pc = color;

    const c = Color.init(@as(u8, @intFromFloat(pc.x() * 255)), @as(u8, @intFromFloat(pc.y() * 255)), @as(u8, @intFromFloat(pc.z() * 255)), @as(u8, @intFromFloat(pc.w() * 255)));

    if (px.x() >= 0 and px.x() <= 1000 and px.y() >= 0 and px.y() <= 1000)
        writePixel(@as(i32, @intFromFloat(px.x())), @as(i32, @intFromFloat(px.y())), px.z(), c);
}

///
pub fn drawWorldLine(mvp: *const Mat44f, start: Vec4f, end: Vec4f, color: Vec4f, shader: *Material) void {
    const spx = shader.vertexShader(mvp, 0, start);
    const epx = shader.vertexShader(mvp, 0, end);
    const pc = color;

    const c = Color.init(@as(u8, @intFromFloat(pc.x() * 255)), @as(u8, @intFromFloat(pc.y() * 255)), @as(u8, @intFromFloat(pc.z() * 255)), @as(u8, @intFromFloat(pc.w() * 255)));

    if (spx.x() >= 0 and spx.x() <= 1000 and spx.y() >= 0 and spx.y() <= 1000 and
        epx.x() >= 0 and epx.x() <= 1000 and epx.y() >= 0 and epx.y() <= 1000)
        drawLine(@as(i32, @intFromFloat(spx.x())), @as(i32, @intFromFloat(spx.y())), @as(i32, @intFromFloat(epx.x())), @as(i32, @intFromFloat(epx.y())), c);
}

pub fn drawProgress(x: i16, y: i16, max_width: f32, value: f32, max_value: f32) void {
    const cs = std.math.clamp(value, 0.0, max_value) / max_value;
    // const cs2 = cs*cs;
    drawLine(x, y, @as(c_int, @intFromFloat(cs * max_width)), y, Color.fromNormal(cs, (1 - cs), 0.2, 1));
}

pub fn writePixelNormal(x: i32, y: i32, z: f32, c: Vec4f) void {
    colorBuffer[currentBuffer].write(x, y, Color.fromNormal(c.x(), c.y(), c.z(), c.w()));
    depthBuffer[currentBuffer].write(x, y, z);
}

pub fn writePixel(x: i32, y: i32, z: f32, c: Color) void {
    colorBuffer[currentBuffer].write(x, y, c);
    depthBuffer[currentBuffer].write(x, y, z);
}

pub fn drawSdfBox(
    proj: *const Mat44f,
    mv: *const Mat44f,
    point: Vec4f,
    color: Vec4f,
    shader: *Material,
    box: Vec4f,
) void {
    //const px = shader.projectionShader(proj, shader.vertexShader(mv, 0, point, shader), viewport, shader);

    //const c = Color.fromNormalVec4f(color);

    const offsets = [_]Vec4f{
        Vec4f.init(0.5, 0.5, 0.5, 0.0),
        Vec4f.init(-0.5, 0.5, 0.5, 0.0),

        Vec4f.init(0.5, 0.5, -0.5, 0.0),
        Vec4f.init(-0.5, 0.5, -0.5, 0.0),

        Vec4f.init(0.5, -0.5, 0.5, 0.0),
        Vec4f.init(-0.5, -0.5, 0.5, 0.0),

        Vec4f.init(0.5, -0.5, -0.5, 0.0),
        Vec4f.init(-0.5, -0.5, -0.5, 0.0),
    };

    var boxMin: Vec4f = Vec4f.init(math.inf(f32), math.inf(f32), math.inf(f32), math.inf(f32));
    var boxMax: Vec4f = Vec4f.init(-math.inf(f32), -math.inf(f32), -math.inf(f32), -math.inf(f32));

    inline for (&offsets) |*offset| {
        const bp = point.addDup(offset.*);
        drawPoint(proj, mv, bp, Vec4f.half(), shader);

        const px = shader.projectionShader(proj, shader.vertexShader(mv, 0, bp, shader), viewport, shader);

        boxMin = Vec4f.min(px, boxMin);
        boxMax = Vec4f.max(px, boxMax);
    }

    const center = shader.vertexShader(mv, 0, point, shader);
    //const size = boxMax.sub(boxMin);

    const half = viewport.scaleDup(0.5);

    // center in viewport
    const boxMaxX = @as(usize, @intFromFloat(@max(0.0, boxMax.x())));
    const boxMaxY = @as(usize, @intFromFloat(@max(0.0, boxMax.y())));

    const boxMinX = @as(usize, @intFromFloat(@max(0.0, boxMin.x())));
    const boxMinY = @as(usize, @intFromFloat(@max(0.0, boxMin.y())));

    _ = box;
    //_ = color;

    for (boxMinY..boxMaxY) |y| {
        for (boxMinX..boxMaxX) |x| {
            const px = Vec4f.init(@as(f32, @floatFromInt(x)), @as(f32, @floatFromInt(y)), boxMax.z() - boxMin.z(), 0);

            // if(depthBuffer[currentBuffer].setLessThan(@intCast(i32,x),@intCast(i32,y),px.z()) == 0)
            //     continue;
            const ray = Vec4f.init(half.x() * px.x() + half.x(), half.y() * -px.y() + half.y(), -1, 0).normalized3();

            const rayDist = ray.dot3(center); //*center.length3();
            const rayPoint = ray.scaleDup(rayDist);
            const boxDist = sdf.Sphere(rayPoint.subDup(center), 0.5);

            if (boxDist >= 0.5)
                continue;

            if (px.x() >= 0 and px.x() <= 1000 and px.y() >= 0 and px.y() <= 1000) {
                writePixel(@as(i32, @intFromFloat(px.x())), @as(i32, @intFromFloat(px.y())), px.z(), Color.fromNormalVec4f(color.scale3Dup(@fabs(boxDist))));
            }
        }
    }
}
