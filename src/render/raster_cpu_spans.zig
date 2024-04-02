const std = @import("std");
const warn = std.debug.print;
const fmt = std.fmt;
const assert = std.debug.assert;
const math = std.math;

const tracy = @import("../tracy.zig");
const trace = tracy.trace;

const core = @import("../core/core.zig");

const Vec4f = core.vector.Vec4f;
const Vec4i = core.vector.Vec4i;
const Mat44f = core.Mat44f;
const Bounds = core.Bounds;
const Color = core.Color;
const Profile = core.Profile;

pub const Font = @import("font.zig").Font;
pub const Mesh = @import("mesh.zig").Mesh;
pub const Material = @import("material.zig").Material;

pub const PixelBuffer = @import("pixel_buffer.zig").PixelBuffer;
pub const PixelRenderer = @import("pixel_buffer.zig").PixelRenderer;
pub const ColorPixelBuffer = PixelBuffer(Color);
pub const DepthPixelBuffer = PixelBuffer(f32);

const RasterData = struct {
    pub const MeshList = std.MultiArrayList(MeshRasterData);
    pub const TriangleList = std.MultiArrayList(TriRasterData);
    pub const TriangleSpanList = std.MultiArrayList(TriSpanData);

    allocator: std.mem.Allocator,

    stats: Stats,
    viewport: Vec4f = Vec4f.zero(),
    renderBounds: Bounds = undefined,

    meshes: MeshList,
    triangles: TriangleList,
    spans: TriangleSpanList,

    colorBuffer: [2]ColorPixelBuffer = undefined,
    depthBuffer: [2]DepthPixelBuffer = undefined,

    currentBuffer: u1 = 0,

    profile: ?*Profile = undefined,
    colorRenderer: PixelRenderer(Color) = undefined,

    pub fn init(_allocator: std.mem.Allocator) RasterData {
        return .{
            .allocator = _allocator,
            .state = Stats.init(),
            .viewport = Vec4f.zero(),
            .renderBounds = Bounds{},
            .meshes = MeshList{},
            .triangles = TriangleList{},
            .spans = TriangleSpanList{},
        };
    }
};

pub const Stats = struct {
    const Self = @This();

    queuedMeshes: u32 = 0,
    totalMeshes: u32 = 0,
    totalTris: u32 = 0,
    totalPixels: u32 = 0,

    renderedMeshes: u32 = 0,
    renderedTris: u32 = 0,
    renderedPixels: usize = 0,

    trisTooSmall: u32 = 0,
    trisTooNear: u32 = 0,
    trisBackfacing: u32 = 0,
    jobWaitCount: u32 = 0,
    jobCount: u32 = 0,

    renderStartTime: u64 = 0,
    renderEndTime: u64 = 0,

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

    // pub fn trace(self: Stats) void {
    //     const ti = @typeInfo(Stats);
    //     inline for (ti.Struct.fields) |field| {
    //         tracy.plotValue(field.name.ptr, @field(self, field.name));
    //     }
    // }

    pub fn print(self: *Self) void {
        std.debug.print("{} ns [ m({}, {}|{d:.2}%), t({}, {}|{d:.2}%, <{d}, |<{d}, bf{d}), p({}, {}|{d:.2}%) ]\n", .{
            (self.renderEndTime - self.renderStartTime),
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

pub const MeshRasterData = struct {
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

pub const TriRasterData = struct {
    pub const Visible = 1 << 0;
    pub const TooSmall = 1 << 1;

    id: u32,
    offset: u16,

    meshData: MeshRasterData = undefined,

    indicies: [3]u16 = undefined,

    color: [3]Vec4f = undefined,
    uv: [3]Vec4f = undefined,
    normals: [3]Vec4f = undefined,
    worldNormals: [3]Vec4f = undefined,

    rawVertex: [3]Vec4f = undefined,

    // transformed and projected verticies
    screenVertex: [3]Vec4f = undefined, // Projection * Model * View
    cameraVertex: [3]Vec4f = undefined, // Model * View
    worldVertex: [3]Vec4f = undefined, // Model

    normal: Vec4f = Vec4f.zero(),
    edges: [3]Vec4f = undefined,
    screenArea: f32 = 0,
    flags: u32 = 0,
    backfacing: bool = false,

    screenBounds: Bounds = undefined,
    ySortedIndex: [3]usize = undefined,

    pixelsFilled: usize = 0,
};

const TriSpanData = struct {
    triData: *TriRasterData,
};

const Plane = struct {
    position: Vec4f = undefined,
    normal: Vec4f = undefined,

    pub fn pointDistance(self: Plane, point: Vec4f) f32 {
        return point.subDup(self.position).dot3(self.normal);
    }
};

const Frustum = struct {
    const Self = @This();

    const Planes = enum(u8) {
        Near = 0,
        Far = 1,
        Right = 2,
        Left = 3,
        Top = 4,
        Bottom = 5,
    };

    planes: [6]Plane = undefined,

    pub fn near(self: Self) Plane {
        return self.planes[@intFromEnum(Planes.Near)];
    }

    pub fn from(self: *Self, view: *const Mat44f, aspect: f32, fovY: f32, znear: f32, zfar: f32) void {
        const PI = 3.1415926535897932384626433832795;
        const DEG2RAD = PI / 180.0;

        const halfV = zfar * @tan(fovY * 0.5 * DEG2RAD);
        const halfH = halfV * aspect;

        const position = view.position();
        const p2 = view.mul_vec4(Vec4f.zero());

        const forward = view.forward();
        const forwardFar = forward.scaleDup(zfar);
        const up = view.up();
        const right = view.right();

        self.planes[@intFromEnum(Planes.Near)] = .{
            .position = position.addDup(forward.scale3Dup(znear)),
            .normal = forward.negDup(),
        };
        self.planes[@intFromEnum(Planes.Far)] = .{
            .position = p2.addDup(forwardFar),
            .normal = forward,
        };

        self.planes[@intFromEnum(Planes.Right)] = .{
            .position = position,
            .normal = forwardFar
                .subDup(right.scale3Dup(halfH)
                .cross3(up)),
        };
        self.planes[@intFromEnum(Planes.Left)] = .{
            .position = position,
            .normal = up.cross3(forwardFar.addDup(right.scale3Dup(halfH))),
        };

        self.planes[@intFromEnum(Planes.Top)] = .{
            .position = position,
            .normal = right.cross3(forwardFar.subDup(up.scale3Dup(halfV))),
        };
        self.planes[@intFromEnum(Planes.Bottom)] = .{
            .position = position,
            .normal = forwardFar.addDup(up.scale3Dup(halfV)).cross3(right),
        };
    }

    pub fn isPointInside(self: Self, point: Vec4f) bool {
        if (self.planes[0].pointDistance(point) < 0.0)
            return false;

        // const len = self.planes.len;
        // inline for (0..len) |p| {
        //     if (self.planes[p].pointDistance(point) < 0.0)
        //         return false;
        // }

        return true;
    }
};

var stats = Stats.init();
var colorBuffer: [2]PixelBuffer(Color) = undefined;
var depthBuffer: [2]PixelBuffer(f32) = undefined;
var currentBuffer: u1 = 0;
var viewport = Vec4f.zero();
var renderBounds: Bounds = undefined;

var profile: ?*Profile = undefined;
var allocator: *std.mem.Allocator = undefined;
var colorRenderer: PixelRenderer(Color) = undefined;
var triQueue: std.MultiArrayList(TriRasterData) = undefined;
pub var viewFrustum: Frustum = Frustum{};

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

    const clearColor = Color.init(34, 34, 34, 255);
    const clearDepth = std.math.inf(f32);

    triQueue = std.MultiArrayList(TriRasterData){};
    try triQueue.ensureTotalCapacity(alloc.*, 1024 * 1024);

    currentBuffer = 0;

    colorBuffer[0] = try PixelBuffer(Color).init(renderWidth, renderHeight, clearColor, allocator);
    depthBuffer[0] = try PixelBuffer(f32).init(renderWidth, renderHeight, clearDepth, allocator);

    colorBuffer[1] = try PixelBuffer(Color).init(renderWidth, renderHeight, clearColor, allocator);
    depthBuffer[1] = try PixelBuffer(f32).init(renderWidth, renderHeight, clearDepth, allocator);

    inline for (&colorBuffer) |*buffer| {
        buffer.clearDefault();
    }

    colorRenderer = PixelRenderer(Color).init(&colorBuffer[currentBuffer]);

    viewport = Vec4f.init(@as(f32, @floatFromInt(colorBuffer[0].w)), @as(f32, @floatFromInt(colorBuffer[0].h)), 0, 0);
    renderBounds = Bounds.init(Vec4f.zero(), viewport);
}

pub fn shutdown() void {
    colorBuffer[0].deinit();
    colorBuffer[1].deinit();
    depthBuffer[0].deinit();
    depthBuffer[1].deinit();
    triQueue.deinit(allocator.*);
}

pub fn beginFrame(profiler: ?*Profile) *u8 {
    const zone = trace(@src());
    defer zone.end();

    profile = profiler;
    var pprof = profile.?.beginSample("render.beginFrame");
    defer profile.?.endSample(pprof);

    stats.reset();
    const now = std.time.Instant.now() catch std.time.Instant{ .timestamp = 0 };
    stats.renderStartTime = now.timestamp;
    currentBuffer = ~currentBuffer;

    colorBuffer[currentBuffer].clearDefault();
    depthBuffer[currentBuffer].clearDefault();

    colorRenderer.setBuffer(&colorBuffer[currentBuffer]);

    const w: u32 = @intCast(colorBuffer[currentBuffer].w);
    const h: u32 = @intCast(colorBuffer[currentBuffer].h);

    stats.totalPixels = w * h;

    return &colorBuffer[currentBuffer].bufferStart().color[0];
}

fn renderTris() void {
    const zone = trace(@src());
    defer zone.end();

    for (0..triQueue.len) |t| {
        var tri = triQueue.get(t);

        //shadeTriSpan(&tri);
        shadeTriBB(&tri);
    }

    triQueue.len = 0;
}

pub fn endFrame() void {
    const now = std.time.Instant.now() catch std.time.Instant{ .timestamp = 0 };
    stats.renderEndTime = now.timestamp;

    renderTris();
    //pixels.swapBuffers();
    {
        const zone = trace(@src());
        defer zone.end();
    }

    stats.print();
}

fn getTriangleNormal(points: [3]Vec4f) Vec4f {
    var e0 = points[1];
    var e1 = points[2];

    e0.sub(points[0]);
    e1.sub(points[0]);

    return e0.cross3(e1).normalized3();
}

pub fn drawMesh(m: *const Mat44f, v: *const Mat44f, p: *const Mat44f, mesh: *Mesh, shader: *Material) !void {
    const zone = trace(@src());
    defer zone.end();

    var data = MeshRasterData{
        .mesh = mesh,
        .model = m,
        .view = v,
        .proj = p,
        .shader = shader,
        .mv = undefined,
        .mvp = undefined,
        .vp = undefined,
        .offset = 0,
    };

    // Model * View
    var mv = v.*;
    mv.mul(m.*);
    data.mv = mv;

    // View * projection
    var vp = p.*;
    vp.mul(v.*);

    // Projection * Model * View
    var mvp = p.*;
    mvp.mul(v.*);
    mvp.mul(m.*);
    data.mvp = mvp;

    var t: u16 = 0;
    const ids = mesh.indexBuffer.len;

    // drawWorldLine(&mvp, Vec4f.zero(), Vec4f.right(), Color.red().toNormalVec4f(), shader);
    // drawWorldLine(&mvp, Vec4f.zero(), Vec4f.up(), Color.green().toNormalVec4f(), shader);
    // drawWorldLine(&mvp, Vec4f.zero(), Vec4f.forward(), Color.blue().toNormalVec4f(), shader);

    _ = @atomicRmw(u32, &stats.totalMeshes, .Add, 1, .SeqCst);

    // renderQueue.pending.lockForWrite();
    // defer renderQueue.pending.unlockForWrite();

    var tri = TriRasterData{
        .id = 0,
        .offset = 0,
        .meshData = data,
        .indicies = undefined,
    };

    var pixelsFilled: usize = 0;
    while (t < ids / 3) {
        tri.pixelsFilled = 0;
        tri.id = @as(u16, @truncate(t));
        tri.offset = t * 3;

        if (drawTriJob(&tri)) {
            triQueue.appendAssumeCapacity(tri);
        }
        t += 1;

        // if (tri.pixelsFilled > 0)
        //     break;
    }

    stats.renderedPixels += pixelsFilled;
    _ = @atomicRmw(u32, &stats.renderedMeshes, .Add, 1, .SeqCst);
}

inline fn applyPixelShader(mv: *const Mat44f, mvp: *const Mat44f, pixel: Vec4f, color: Vec4f, normal: Vec4f, uv: Vec4f, material: *Material) Vec4f {
    const zone2 = trace(@src());
    defer zone2.end();
    // var c = color.addDup(
    //   engine.Vec4f.init(
    //     (std.math.sin(uv.x*uv.y*1000)+1/2),
    //     (std.math.cos(uv.y*1000)+1/2),
    //     0,1)
    //   );

    // _ = color;
    _ = pixel;
    _ = mvp;
    // _ = uv;
    _ = mv;

    var c = material.texture.sampleBilinear(uv.x(), uv.y());

    c.lerp(color, 0.5);
    c.setW(1.0);
    //var c = material.texture.sampleBilinear(uv.x(), uv.y());

    //var c = color; //material.texture.sample(uv.x(), uv.y());

    const l = @max(normal.dot3(material.lightDirection) * material.lightIntensity, 0.4);

    c.scale(l);

    return c;
    //return uncharted2_filmic(c);
    //kreturn reinhard(c);
    //return c.scaleDup(l);
}

fn drawTriJob(triJob: *TriRasterData) bool {
    const zone = trace(@src());
    defer zone.end();

    var tri = triJob;
    const data = tri.meshData;
    const mesh = data.mesh;
    const shader = data.shader;

    _ = @atomicRmw(u32, &stats.totalTris, .Add, 1, .SeqCst);

    comptime var p = 0;
    var insideCount: i32 = 0;
    inline while (p < 3) {
        var offset = tri.offset + p;
        assert(offset < mesh.indexBuffer.len);
        tri.indicies[p] = mesh.indexBuffer[offset];
        tri.rawVertex[p] = mesh.vertexBuffer[tri.indicies[p]];

        tri.cameraVertex[p] = shader.vertexShader(&data.mv, offset, tri.rawVertex[p], shader);
        tri.screenVertex[p] = shader.projectionShader(data.proj, tri.cameraVertex[p], viewport, shader);

        tri.normals[p] = mesh.vertexNormalBuffer[mesh.indexNormalBuffer[offset]];
        tri.worldNormals[p] = data.model.mul33_vec4(tri.normals[p]);
        tri.uv[p] = Vec4f.init(mesh.textureCoordBuffer[mesh.indexUVBuffer[(offset)] * 2 + 0], mesh.textureCoordBuffer[mesh.indexUVBuffer[(offset)] * 2 + 1], 0, 0);

        const w = tri.screenVertex[p].w();
        if (w != 0.0) {
            tri.uv[p].div(w);
            //tri.screenVertex[p].divVec(Vec4f.init(1.0, 1.0, w, 1.0));
            tri.screenVertex[p].div(w);
            tri.uv[p].setW(1.0 / w);
        }

        const inside = tri.screenVertex[p].abs().maxElement() <= 1.0;
        insideCount += @intFromBool(inside);

        const half = viewport.scaleDup(0.5);
        const out = tri.screenVertex[p];
        tri.screenVertex[p].setX(half.x() * out.x() + half.x());
        tri.screenVertex[p].setY(half.y() * -out.y() + half.y());

        // tri.screenVertex[p].addScalar(1.0);
        // tri.screenVertex[p].div(2);

        tri.worldVertex[p] = data.model.mul_vec4(tri.rawVertex[p]);

        // const inside = @intFromBool(viewFrustum.isPointInside(tri.worldVertex[p]));

        //const indexScalar = @as(f32, @floatFromInt(offset)) / @as(f32, @floatFromInt(mesh.indexBuffer.len));
        //tri.color[p] = Vec4f.init(0.4 + @as(f32, @floatFromInt(1 - inside)), 0.7, 0.5, 1.0).mulDup(Vec4f.init(indexScalar, indexScalar, indexScalar, 1.0)); //mesh.colorBuffer[tri.indicies[p]];

        tri.color[p] = mesh.colorBuffer[tri.indicies[p]];

        // tri.color[p].setX(tri.color[p].x() + @as(f32, @floatFromInt(1 - inside)));
        p += 1;
    }

    // frustum cull
    if (insideCount < 1) {
        return false;
    }

    tri.normal = getTriangleNormal(tri.cameraVertex);
    tri.backfacing = shader.backfaceCull and tri.cameraVertex[0].normalized3().dot3(tri.normal) > 0.0000000000001;
    tri.screenArea = Vec4f.triArea(tri.screenVertex[0], tri.screenVertex[1], tri.screenVertex[2]);
    tri.screenBounds = Bounds.init(tri.screenVertex[0], tri.screenVertex[0]);

    if (tri.backfacing) {
        _ = @atomicRmw(u32, &stats.trisBackfacing, .Add, 1, .SeqCst);
        return false;
    }

    p = 1;
    inline while (p < 3) {
        tri.screenBounds.add(tri.screenVertex[p]);
        p += 1;
    }

    if (std.math.isInf(tri.screenBounds.min.maxElement()) or std.math.isInf(tri.screenBounds.max.maxElement())) {
        return false;
    }

    // // Too small to see
    // if (tri.screenArea <= 0) {
    //     _ = @atomicRmw(u32, &stats.trisTooSmall, .Add, 1, .SeqCst);
    //     return false;
    // }

    const rb = renderBounds;
    tri.screenBounds.limit(rb);
    tri.screenBounds.topLeftHandLimit();

    tri.ySortedIndex = .{ 0, 1, 2 };

    var t0 = tri.screenVertex[0];
    var t1 = tri.screenVertex[1];
    var t2 = tri.screenVertex[2];

    // sort by y value
    if (t0.y() > t1.y()) {
        std.mem.swap(Vec4f, &t0, &t1);
        std.mem.swap(usize, &tri.ySortedIndex[0], &tri.ySortedIndex[1]);
    }
    if (t0.y() > t2.y()) {
        std.mem.swap(Vec4f, &t0, &t2);
        std.mem.swap(usize, &tri.ySortedIndex[0], &tri.ySortedIndex[2]);
    }
    if (t1.y() > t2.y()) {
        std.mem.swap(Vec4f, &t1, &t2);
        std.mem.swap(usize, &tri.ySortedIndex[1], &tri.ySortedIndex[2]);
    }

    // drawWorldLine(&data.mvp, tri.worldVertex[0], tri.worldVertex[1], tri.color[0], shader);
    // drawWorldLine(&data.mvp, tri.worldVertex[1], tri.worldVertex[2], tri.color[1], shader);
    // drawWorldLine(&data.mvp, tri.worldVertex[2], tri.worldVertex[0], tri.color[2], shader);

    // drawLine(@as(i32, @intFromFloat(tri.screenVertex[0].x())), @as(i32, @intFromFloat(tri.screenVertex[0].y())), @as(i32, @intFromFloat(tri.screenVertex[1].x())), @as(i32, @intFromFloat(tri.screenVertex[1].y())), Color.fromNormalVec4f(tri.color[0]));

    // drawLine(@as(i32, @intFromFloat(tri.screenVertex[1].x())), @as(i32, @intFromFloat(tri.screenVertex[1].y())), @as(i32, @intFromFloat(tri.screenVertex[2].x())), @as(i32, @intFromFloat(tri.screenVertex[2].y())), Color.fromNormalVec4f(tri.color[1]));

    // drawLine(@as(i32, @intFromFloat(tri.screenVertex[2].x())), @as(i32, @intFromFloat(tri.screenVertex[2].y())), @as(i32, @intFromFloat(tri.screenVertex[0].x())), @as(i32, @intFromFloat(tri.screenVertex[0].y())), Color.fromNormalVec4f(tri.color[2]));

    // // math.max(math.max(v[0].y(), v[1].y()), v[2].y());

    _ = @atomicRmw(u32, &stats.renderedTris, .Add, 1, .SeqCst);

    return true;
}

fn shadeTriSpan(triData: *TriRasterData) void {
    // const zone = trace(@src());
    // defer zone.end();

    var t0 = triData.screenVertex[triData.ySortedIndex[0]];
    var t1 = triData.screenVertex[triData.ySortedIndex[1]];
    var t2 = triData.screenVertex[triData.ySortedIndex[2]];

    const width: f32 = @floatFromInt(colorBuffer[currentBuffer].w);
    const height: f32 = @floatFromInt(colorBuffer[currentBuffer].h);

    var totalHeight = std.math.clamp(t2.y() - t0.y(), 0, height);

    if (totalHeight <= 0)
        return;

    var uv0 = triData.uv[triData.ySortedIndex[0]];
    var uv1 = triData.uv[triData.ySortedIndex[1]];
    var uv2 = triData.uv[triData.ySortedIndex[2]];

    var vn0 = triData.worldNormals[triData.ySortedIndex[0]];
    var vn1 = triData.worldNormals[triData.ySortedIndex[1]];
    var vn2 = triData.worldNormals[triData.ySortedIndex[2]];

    var vc0 = triData.color[triData.ySortedIndex[0]];
    var vc1 = triData.color[triData.ySortedIndex[1]];
    var vc2 = triData.color[triData.ySortedIndex[2]];

    const dy1 = @fabs(t1.y() - t0.y());
    const dy2 = @fabs(t2.y() - t1.y());
    const flatTop = t1.y() == t0.y();
    const shader = triData.meshData.shader;

    var x: f32 = t0.x();
    var y: f32 = t0.y();
    var z: f32 = t0.z();

    while (y <= t2.y()) {
        defer y += 1;

        const zone2 = trace(@src());
        defer zone2.end();

        const starty = y - t0.y();

        const second_half = starty > dy1 or flatTop;
        const segment_height = @fabs(if (second_half) dy2 else dy1);

        const alpha = starty / totalHeight;
        const beta = (starty - (if (second_half) dy1 else 0)) / segment_height;

        // screen position
        var A = t0.addDup(t2.subDup(t0).scaleDup(alpha));
        var B = if (second_half)
            t1.addDup(t2.subDup(t1).scaleDup(beta))
        else
            t0.addDup(t1.subDup(t0).scaleDup(beta));

        // texture coords
        var U = uv0.addDup(uv2.subDup(uv0).scaleDup(alpha));
        var V = if (second_half)
            uv1.addDup(uv2.subDup(uv1).scaleDup(beta))
        else
            uv0.addDup(uv1.subDup(uv0).scaleDup(beta));

        // vertex normals
        var VA = vn0.addDup(vn2.subDup(vn0).scaleDup(alpha));
        var VB = if (second_half)
            vn1.addDup(vn2.subDup(vn1).scaleDup(beta))
        else
            vn0.addDup(vn1.subDup(vn0).scaleDup(beta));

        // vertex colors
        var CA = vc0.addDup(vc2.subDup(vc0).scaleDup(alpha));
        var CB = if (second_half)
            vc1.addDup(vc2.subDup(vc1).scaleDup(beta))
        else
            vc0.addDup(vc1.subDup(vc0).scaleDup(beta));

        A.clampVec(
            Vec4f.init(0.0, 0.0, A.z(), A.w()),
            Vec4f.init(width - 1, height - 1, A.z(), A.w()),
        );
        B.clampVec(
            Vec4f.init(0.0, 0.0, B.z(), B.w()),
            Vec4f.init(width - 1, height - 1, B.z(), B.w()),
        );

        if (@fabs((A.x())) > @fabs((B.x()))) {
            std.mem.swap(Vec4f, &A, &B);
            std.mem.swap(Vec4f, &U, &V);
            std.mem.swap(Vec4f, &VA, &VB);
            std.mem.swap(Vec4f, &CA, &CB);
        }

        const ix: i32 = @intFromFloat(std.math.clamp(x, 0, width));
        const iy: i32 = @intFromFloat(std.math.clamp(y, 0, height));

        const ldiff = B.subDup(A);

        const tstep: f32 = 1.0 / ldiff.x();
        var t: f32 = 0;
        const depthTest = shader.depthTest;

        var dsx = depthBuffer[currentBuffer].pxIndex(ix, iy);
        var csx = colorBuffer[currentBuffer].pxIndex(ix, iy);

        x = A.x();
        while (x <= B.x()) {
            // const pixtracy = trace(@src());
            // defer pixtracy.end();
            // var pprof = profile.?.beginSample("render.mesh.draw.tri.pixel");
            // defer profile.?.endSample(pprof);
            //_ = @atomicRmw(u32, &stats.totalPixels, .Add, 1, .SeqCst);
            defer x += 1;

            const lerpZone = tracy.traceNamed(@src(), "lerpZone");

            //const x = @as(i32, @intCast(j));
            const uvs = U.lerpDup(V, t);
            const pcuvs = uvs.divDup(uvs.w());
            const vns = VA.lerpDup(VB, t);
            const vc = CA.lerpDup(CB, t);
            const p = A.lerpDup(B, t);

            z = p.z();
            lerpZone.end();

            var write = true;
            if (depthTest == 1) {
                const depthZone = tracy.traceNamed(@src(), "depthTest");
                defer depthZone.end();
                const depth = depthBuffer[currentBuffer].readIndex(dsx);
                write = std.math.isInf(depth) or depth < z;
            }

            if (write) {
                const c = shader.pixelShader(
                    &triData.meshData.mv,
                    &triData.meshData.mvp,
                    p,
                    vc,
                    vns,
                    pcuvs,
                    shader,
                );
                // const c = applyPixelShader(&triData.meshData.mv, &triData.meshData.mvp, p, vc, vns, pcuvs, shader);

                triData.pixelsFilled += 1;
                //writePixel(x, y, z, Color.fromNormalVec4f(c));

                colorBuffer[currentBuffer].writeIndex(csx, Color.fromNormalVec4f(c));
                depthBuffer[currentBuffer].writeIndex(dsx, z);
            }
            t += tstep;
            dsx += 1;
            csx += 1;
        }
    }
}

fn shadeTriBB(triData: *TriRasterData) void {
    const bounds = triData.screenBounds;
    var y = bounds.min.y();
    var p: Vec4f = Vec4f.init(0, 0, 0, 0);
    var pixelNormal: Vec4f = Vec4f.init(0, 0, 0, 0);
    var fbc: Vec4f = Vec4f.init(0, 0, 0, 1);
    var uv: Vec4f = Vec4f.init(0, 0, 0, 1);
    // var c: Color = Color.black();
    const shader = triData.meshData.shader;

    _ = @atomicRmw(u32, &stats.renderedTris, .Add, 1, .SeqCst);

    const minx: i32 = @intFromFloat(bounds.min.x());
    const miny: i32 = @intFromFloat(bounds.min.y());

    var dsy = depthBuffer[currentBuffer].pxIndex(minx, miny);
    const width: usize = @intCast(colorBuffer[currentBuffer].w);
    const depthTest = shader.depthTest;

    while (y <= bounds.max.y()) {
        var x = bounds.min.x();
        defer y += 1;
        defer dsy += width;

        var dsx = dsy;
        while (x <= bounds.max.x()) {
            // const pixtracy = trace(@src());
            // defer pixtracy.end();
            // var pprof = profile.?.beginSample("render.mesh.draw.tri.pixel");
            // defer profile.?.endSample(pprof);
            //_ = @atomicRmw(u32, &stats.totalPixels, .Add, 1, .SeqCst);
            defer dsx += 1;
            defer x += 1;

            p.setX(x);
            p.setY(y);

            var triBary = Vec4f.triBarycentericCoordsOld(triData.screenVertex[0], triData.screenVertex[1], triData.screenVertex[2], p);

            if (std.math.signbit(triBary.minElement()))
                continue;

            triBary.div(triData.screenArea);

            const z = (triBary.x() * triData.screenVertex[0].z() +
                triBary.y() * triData.screenVertex[1].z() +
                triBary.z() * triData.screenVertex[2].z());

            p.setZ(z);

            const pt = triBary.triInterpArray(triData.screenVertex, 1.0, 1.0);

            var write = true;
            if (depthTest == 1) {
                // const depthZone = tracy.traceNamed(@src(), "depthTest");
                // defer depthZone.end();
                const depth = depthBuffer[currentBuffer].readIndex(dsx);
                write = std.math.isInf(depth) or depth > z;
            }

            if (write) {
                // interpolate vertex colors across all pixels
                fbc = triBary.triInterpArray(triData.color, 1.0, 1.0);
                pixelNormal = triBary.triInterpArray(triData.worldNormals, 1.0, 1.0);
                uv = triBary.triInterpArray(triData.uv, 1.0, 1.0);
                uv.div(uv.w());

                const vc = shader.pixelShader(&triData.meshData.mv, &triData.meshData.mvp, pt, fbc, pixelNormal, uv, shader);

                triData.pixelsFilled += 1;

                colorBuffer[currentBuffer].writeIndex(dsx, Color.fromNormalVec4f(vc));
                depthBuffer[currentBuffer].writeIndex(dsx, z);
            }
        }
    }
}

pub fn drawPointMesh(mvp: *const Mat44f, mesh: *Mesh, shader: *Material) void {
    // const ids = mesh.vertexBuffer.len;
    _ = shader;
    for (mesh.vertexBuffer, 0..) |vertex, i| {
        drawPoint(mvp, vertex, mesh.colorBuffer[i]);
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

        var coy: i32 = 0;
        while (coy < font.glyphHeight) {
            defer coy += 1;

            var cox: i32 = 0;
            while (cox < font.glyphWidth) {
                defer cox += 1;

                var samp = font.characterColor(cx, cy, cox, coy);
                if (samp.x() <= 0.001)
                    continue;

                writePixelNormal((x + (font.glyphWidth - 1) * @as(i32, @intCast(drawCount))) + cox, starty + coy, 1.0, samp);
            }
        }

        drawCount += 1;
    }
}

//
pub fn drawPoint(mvp: *const Mat44f, point: Vec4f, color: Vec4f, shader: *Material) void {
    const px = shader.vertexShader(mvp, 0, point, shader);
    const pc = color;

    const c = Color.init(@as(u8, @intFromFloat(pc.x() * 255)), @as(u8, @intFromFloat(pc.y() * 255)), @as(u8, @intFromFloat(pc.z() * 255)), @as(u8, @intFromFloat(pc.w() * 255)));

    if (px.x() >= 0 and px.x() <= 1000 and px.y() >= 0 and px.y() <= 1000)
        writePixel(@as(i32, @intFromFloat(px.x())), @as(i32, @intFromFloat(px.y())), c);
}

///
pub fn drawWorldLine(mvp: *const Mat44f, start: Vec4f, end: Vec4f, color: Vec4f, shader: *Material) void {
    const spx = shader.vertexShader(mvp, 0, start, shader);
    const epx = shader.vertexShader(mvp, 0, end, shader);
    // const pc = color;

    _ = color;
    const c = Color.white(); //Color.init(@as(u8, @intFromFloat(pc.x() * 255)), @as(u8, @intFromFloat(pc.y() * 255)), @as(u8, @intFromFloat(pc.z() * 255)), @as(u8, @intFromFloat(pc.w() * 255)));
    const screen = renderBounds.size();

    if (spx.x() >= 0 and spx.x() < screen.x() and spx.y() >= 0 and spx.y() <= screen.x() and
        epx.x() >= 0 and epx.x() < screen.x() and epx.y() >= 0 and epx.y() <= screen.x())
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
