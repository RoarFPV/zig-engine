const std = @import("std");
const warn = std.debug.print;
const fmt = std.fmt;
const assert = std.debug.assert;
const math = std.math;

const tracy = @import("../tracy.zig");
const trace = tracy.trace;

const core = @import("../core/core.zig");
const interp = @import("../core/interp.zig");

const Vec4f = core.vector.Vec4f;
const Vec4i = core.vector.Vec4i;
const Mat44f = core.Mat44f;
const Bounds = core.Bounds;
const Color = core.Color;
const Profile = core.Profile;

pub const Font = @import("font.zig").Font;
pub const Mesh = @import("mesh.zig").Mesh;
pub const Material = @import("material.zig").Material;
pub const Terrain = @import("terrain.zig").Terrain;

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

    renderStartTime: std.time.Instant = undefined,
    renderEndTime: std.time.Instant = undefined,

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
            @field(self, field.name) = undefined;
        }
    }

    // pub fn trace(self: Stats) void {
    //     const ti = @typeInfo(Stats);
    //     inline for (ti.Struct.fields) |field| {
    //         tracy.plotValue(field.name.ptr, @field(self, field.name));
    //     }
    // }

    pub fn print(self: *Self) void {
        std.debug.print("{d:.3} ms [ m({}, {}|{d:.2}%), t({}, {}|{d:.2}%, <{d}, |<{d}, bf{d}), p({}, {}|{d:.2}%) ]\n", .{
            @as(f32, @floatFromInt(self.renderEndTime.since(self.renderStartTime))) / 1_000_000.0,
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
    const Self = @This();
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
    clipVertex: [3]Vec4f = undefined, // Projection * Model * View
    screenVertex: [3]Vec4f = undefined, // Projection * Model * View / w
    screenTriBary: [4]Vec4f = undefined, // Projection * Model * View / w
    cameraVertex: [3]Vec4f = undefined, // Model * View
    worldVertex: [3]Vec4f = undefined, // Model

    normal: Vec4f = Vec4f.zero(),
    edges: [3]Vec4f = undefined,
    screenArea: f32 = 0,
    flags: u32 = 0,
    backfacing: bool = false,

    screenBounds: Bounds = undefined,
    // ySortedIndex: [3]usize = undefined,

    pixelsFilled: usize = 0,

    clipCount: usize = 0,

    pub fn calcScreenBounds(self: *Self) void {
        self.screenArea = Vec4f.triArea(self.screenVertex[0], self.screenVertex[1], self.screenVertex[2]);
        self.screenBounds = Bounds.init(self.screenVertex[0], self.screenVertex[0]);

        inline for (1..3) |tp| {
            self.screenBounds.add(self.screenVertex[tp]);
        }
    }

    pub fn limitSreenBounds(self: *Self, other: Bounds) void {
        self.screenBounds.limit(other);
        self.screenBounds.topLeftHandLimit();
    }
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
            .normal = forward.neg(),
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

pub var viewConfig: Vec4f = Vec4f.init(0.1, 1000, 90, 0);

pub fn frameStats() Stats {
    return stats;
}

pub fn drawLine(xFrom: i32, yFrom: i32, xTo: i32, yTo: i32, color: Color) void {
    const zone = trace(@src());
    defer zone.end();
    colorRenderer.drawLine(xFrom, yFrom, xTo, yTo, color);
}

pub fn drawLineTri(tri: *TriRasterData, color: Color) void {
    const zone = trace(@src());
    defer zone.end();

    const t0 = core.vector.Vec4fToVec4i(tri.screenVertex[0]);
    const t1 = core.vector.Vec4fToVec4i(tri.screenVertex[1]);
    const t2 = core.vector.Vec4fToVec4i(tri.screenVertex[2]);

    colorRenderer.drawLine(t0.x(), t0.y(), t1.x(), t1.y(), color);
    colorRenderer.drawLine(t1.x(), t1.y(), t2.x(), t2.y(), color);
    colorRenderer.drawLine(t2.x(), t2.y(), t0.x(), t0.y(), color);
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

    viewport = Vec4f.init(
        @as(f32, @floatFromInt(colorBuffer[0].w)),
        @as(f32, @floatFromInt(colorBuffer[0].h)),
        0,
        0,
    );
    renderBounds = Bounds.init(
        Vec4f.zero(),
        viewport.subDup(Vec4f.init(1, 1, 0, 0)),
    );
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
    const pprof = profile.?.beginSample("render.beginFrame", 0);
    defer profile.?.endSample(pprof);

    stats.reset();
    const now = std.time.Instant.now() catch std.time.Instant{ .timestamp = 0 };
    stats.renderStartTime = now;
    currentBuffer = ~currentBuffer;

    colorBuffer[currentBuffer].clearDefault();
    depthBuffer[currentBuffer].clearDefault();

    colorRenderer.setBuffer(&colorBuffer[currentBuffer]);

    const w: u32 = @intCast(colorBuffer[currentBuffer].w);
    const h: u32 = @intCast(colorBuffer[currentBuffer].h);

    stats.totalPixels = w * h;

    return &colorBuffer[currentBuffer].bufferStart().color[0];
}

pub var useSingleBB: bool = true;

fn renderTris() void {
    const zone = trace(@src());
    defer zone.end();

    var next: usize = 0;
    while (next < triQueue.len) {
        defer next += 1;

        var tri = triQueue.get(next);

        if (!clipTri(&tri)) {
            continue;
        }
        //shadeTriSpan(&tri);
        // shadeTriBB(&tri, tri.screenBounds);
        shadeTriBB(&tri, tri.screenBounds);
    }

    triQueue.len = 0;
}

fn renderTrisSingle() void {
    const zone = trace(@src());
    defer zone.end();

    var next: usize = 0;
    while (next < triQueue.len) {
        defer next += 1;

        var tri = triQueue.get(next);

        if (!clipTri(&tri)) {
            continue;
        }

        shadeTriBBSingle(&tri, tri.screenBounds);
    }

    triQueue.len = 0;
}

pub fn endFrame() void {
    const now = std.time.Instant.now() catch std.time.Instant{ .timestamp = 0 };
    stats.renderEndTime = now;

    if (useSingleBB) {
        renderTrisSingle();
    } else {
        renderTris();
    }
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

    _ = @atomicRmw(u32, &stats.totalMeshes, .Add, 1, .seq_cst);

    // renderQueue.pending.lockForWrite();
    // defer renderQueue.pending.unlockForWrite();

    var tri = TriRasterData{
        .id = 0,
        .offset = 0,
        .meshData = data,
        .indicies = undefined,
    };

    // var pixelsFilled: usize = 0;
    while (t < ids / 3) {
        defer t += 1;
        tri.pixelsFilled = 0;
        tri.id = @as(u16, @truncate(t));
        tri.offset = t * 3;

        if (drawTriJob(&tri)) {
            triQueue.appendAssumeCapacity(tri);
        }

        // if (tri.pixelsFilled > 0)
        // break;
    }

    // stats.renderedPixels += pixelsFilled;
    _ = @atomicRmw(u32, &stats.renderedMeshes, .Add, 1, .seq_cst);
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

const IntersectClipResult = struct {
    v: Vec4f,
    delta: f32,
};

fn intersectClipEdge(outPoint: Vec4f, inPoint: Vec4f, axis: usize) IntersectClipResult {
    const w = outPoint.w();
    // const aval = p.v[axis];
    const dw = @abs(w - outPoint.v[axis]);

    const toOut = outPoint.subDup(inPoint);

    const total = toOut.v[axis];
    const delta = 1 - dw / total;

    // point along this edge that lies on boundary
    const offset = toOut.scaleDup(delta);

    var np = inPoint.addDup(offset);
    np.v[axis] = w;

    return .{ .v = np, .delta = delta };
}

fn addClippedTri(tri: *const TriRasterData, v0: Vec4f, v1: Vec4f, v2: Vec4f) void {
    var tri1 = tri.*;
    tri1.clipCount += 1;

    tri1.cameraVertex[0] = v0;
    tri1.cameraVertex[1] = v1;
    tri1.cameraVertex[2] = v2;

    const shader = tri1.meshData.shader;
    inline for (0..3) |i| {
        tri1.cameraVertex[i].setW(1);
        tri1.clipVertex[i] = shader.projectionShader(
            tri1.meshData.proj,
            tri1.cameraVertex[i],
            viewport,
            shader,
        );

        const clipV = tri1.clipVertex[i];
        tri1.screenVertex[i] = clipV.divDup(clipV.w());

        const half = viewport.scaleDup(0.5);
        const out = tri1.screenVertex[i];
        tri1.screenVertex[i].setX(half.x() * out.x() + half.x());
        tri1.screenVertex[i].setY(half.y() * -out.y() + half.y());
    }

    tri1.calcScreenBounds();
    tri1.limitSreenBounds(renderBounds);

    // TODO: interpolate all values
    triQueue.appendAssumeCapacity(tri1);
}

fn clipTri(tri: *TriRasterData) bool {
    var srt = core.profiler.Sampler.initAndBegin(profile.?, @src().fn_name, 2);
    defer srt.end();

    var totalOutside: usize = 0;

    if (tri.clipCount > 0)
        return true;

    var clipped = [3]Vec4f{
        tri.cameraVertex[0], //.divDup(tri.clipVertex[0].w()),
        tri.cameraVertex[1], //.divDup(tri.clipVertex[1].w()),
        tri.cameraVertex[2], //.divDup(tri.clipVertex[2].w()),
    };

    const nearZ = -viewConfig.x();
    // inline for (0..3) |axis| {
    const axis = 2; // z only
    inline for (0..3) |v| {
        const inside = clipped[v].v[axis] <= nearZ;
        tri.edges[v].setW(@floatFromInt(@intFromBool(inside)));
        totalOutside += @intFromBool(!inside);
    }

    if (totalOutside == 0)
        return true;

    // all points outside of a single plane?
    if (totalOutside == 3)
        return false;

    inline for (0..3) |v| {
        clipped[v].setW(nearZ);
    }

    const nextIndex = .{ 1, 2, 0 };
    const prevIndex = .{ 2, 0, 1 };

    inline for (0..3) |v| {
        const next = nextIndex[v];
        const prev = prevIndex[v];

        const p = clipped[v];

        if (tri.edges[v].w() > 0) {
            // point along this edge that lies on boundary

            const np = clipped[next];
            const pp = clipped[prev];

            const nextOutside = (np.v[axis]) > nearZ;
            const prevOutside = (pp.v[axis]) > nearZ;

            const npi = intersectClipEdge(np, p, axis);
            const ppi = intersectClipEdge(pp, p, axis);

            if (nextOutside and prevOutside) {
                // both out side, just need to make one tri
                // TODO: need to check winding order
                addClippedTri(tri, npi.v, ppi.v, p);
                return false;
            }

            // only one outside
            // need to make two new tris

            if (nextOutside) {
                const out2pi = intersectClipEdge(np, pp, axis);

                addClippedTri(tri, p, npi.v, pp);
                addClippedTri(tri, npi.v, out2pi.v, pp);
            } else {
                const out2pi = intersectClipEdge(pp, np, axis);

                addClippedTri(tri, p, np, ppi.v);
                addClippedTri(tri, ppi.v, np, out2pi.v);
            }

            return false;
        }
    }
    // }

    return true;
}

fn drawTriJob(triJob: *TriRasterData) bool {
    const zone = trace(@src());
    defer zone.end();

    var srt = core.profiler.Sampler.initAndBegin(profile.?, @src().fn_name, 1);
    defer srt.end();

    var tri = triJob;
    const data = tri.meshData;
    const mesh = data.mesh;
    const shader = data.shader;

    _ = @atomicRmw(u32, &stats.totalTris, .Add, 1, .seq_cst);

    inline for (0..3) |p| {
        var offset: u16 = @intCast(p);
        offset += tri.offset;

        assert(offset < mesh.indexBuffer.len);
        tri.indicies[p] = mesh.indexBuffer[offset];
        tri.rawVertex[p] = mesh.vertexBuffer[tri.indicies[p]];

        tri.cameraVertex[p] = shader.vertexShader(&data.mv, offset, tri.rawVertex[p], shader);
        tri.clipVertex[p] = shader.projectionShader(data.proj, tri.cameraVertex[p], viewport, shader);

        tri.normals[p] = mesh.vertexNormalBuffer[mesh.indexNormalBuffer[offset]];
        tri.worldNormals[p] = data.model.mul33_vec4(tri.normals[p]).normalized3();
        tri.uv[p] = Vec4f.init(
            mesh.textureCoordBuffer[mesh.indexUVBuffer[(offset)] * 2 + 0],
            mesh.textureCoordBuffer[mesh.indexUVBuffer[(offset)] * 2 + 1],
            0,
            0,
        );

        const n = tri.clipVertex[p].maxElement();
        if (std.math.isInf(n) or std.math.isNan(n)) {
            var a = tri.clipVertex[p];
            _ = a.addDup(a);
        }

        const w = tri.clipVertex[p].w();

        // tri.clipVertex[p].abs().

        // const inside = tri.clipVertex[p].abs().maxElement() <= std.math.fabs(w);

        // insideCount += @intFromBool(inside);

        if (w != 0.0) {
            tri.uv[p].div(w);
            tri.worldNormals[p].div(w);
            tri.color[p].div(w);
            tri.screenVertex[p] = tri.clipVertex[p].divDup(w);
            tri.edges[p] = tri.clipVertex[p].divDup(w);
            tri.uv[p].setW(1 / w);
        }

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
    }

    // frustum cull
    // if (insideCount == 0) {
    //     return false;
    // }
    tri.normal = getTriangleNormal(tri.cameraVertex);
    tri.backfacing = shader.backfaceCull and tri.cameraVertex[0].normalized3().dot3(tri.normal) > 0.0000000000001;

    if (tri.backfacing) {
        _ = @atomicRmw(u32, &stats.trisBackfacing, .Add, 1, .seq_cst);
        return false;
    }

    // inline for (0..3) |axis| {
    //     var insideCount: usize = 0;
    //     inline for (0..3) |v| {
    //         insideCount += @intFromBool(std.math.fabs(tri.clipVertex[v].v[axis]) <= tri.clipVertex[v].w());
    //     }

    //     // frustum culling
    //     if (insideCount == 0)
    //         return false;
    // }

    tri.edges[0] = (tri.clipVertex[1].subDup(tri.clipVertex[0]));
    tri.edges[1] = (tri.clipVertex[2].subDup(tri.clipVertex[1]));
    tri.edges[2] = (tri.clipVertex[0].subDup(tri.clipVertex[2]));

    inline for (0..4) |i| {
        tri.screenTriBary[i] = Vec4f.init(
            tri.screenVertex[0].v[i],
            tri.screenVertex[1].v[i],
            tri.screenVertex[2].v[i],
            0,
        );
    }

    tri.calcScreenBounds();

    // if (std.math.isInf(tri.screenBounds.min.maxElement()) or std.math.isInf(tri.screenBounds.max.maxElement())) {
    //     return false;
    // }

    // // Too small to see
    if (@abs(tri.screenArea) <= 0) {
        _ = @atomicRmw(u32, &stats.trisTooSmall, .Add, 1, .seq_cst);
        return false;
    }

    tri.limitSreenBounds(renderBounds);

    // tri.ySortedIndex = .{ 0, 1, 2 };

    // var t0 = tri.screenVertex[0];
    // var t1 = tri.screenVertex[1];
    // var t2 = tri.screenVertex[2];

    // // sort by y value
    // if (t0.y() > t1.y()) {
    //     std.mem.swap(Vec4f, &t0, &t1);
    //     std.mem.swap(usize, &tri.ySortedIndex[0], &tri.ySortedIndex[1]);
    // }
    // if (t0.y() > t2.y()) {
    //     std.mem.swap(Vec4f, &t0, &t2);
    //     std.mem.swap(usize, &tri.ySortedIndex[0], &tri.ySortedIndex[2]);
    // }
    // if (t1.y() > t2.y()) {
    //     std.mem.swap(Vec4f, &t1, &t2);
    //     std.mem.swap(usize, &tri.ySortedIndex[1], &tri.ySortedIndex[2]);
    // }

    // drawWorldLine(&data.mvp, tri.worldVertex[0], tri.worldVertex[1], tri.color[0], shader);
    // drawWorldLine(&data.mvp, tri.worldVertex[1], tri.worldVertex[2], tri.color[1], shader);
    // drawWorldLine(&data.mvp, tri.worldVertex[2], tri.worldVertex[0], tri.color[2], shader);

    // drawLine(@as(i32, @intFromFloat(tri.screenVertex[0].x())), @as(i32, @intFromFloat(tri.screenVertex[0].y())), @as(i32, @intFromFloat(tri.screenVertex[1].x())), @as(i32, @intFromFloat(tri.screenVertex[1].y())), Color.fromNormalVec4f(tri.color[0]));

    // drawLine(@as(i32, @intFromFloat(tri.screenVertex[1].x())), @as(i32, @intFromFloat(tri.screenVertex[1].y())), @as(i32, @intFromFloat(tri.screenVertex[2].x())), @as(i32, @intFromFloat(tri.screenVertex[2].y())), Color.fromNormalVec4f(tri.color[1]));

    // drawLine(@as(i32, @intFromFloat(tri.screenVertex[2].x())), @as(i32, @intFromFloat(tri.screenVertex[2].y())), @as(i32, @intFromFloat(tri.screenVertex[0].x())), @as(i32, @intFromFloat(tri.screenVertex[0].y())), Color.fromNormalVec4f(tri.color[2]));

    // // math.max(math.max(v[0].y(), v[1].y()), v[2].y());

    _ = @atomicRmw(u32, &stats.renderedTris, .Add, 1, .seq_cst);

    return true;
}

fn shadeTriBB(triData: *TriRasterData, bounds: Bounds) void {
    const size = bounds.size();
    const b = bounds;
    // b.add(bounds.min.subDup(Vec4f.one()));
    // b.add(bounds.max.addDup(Vec4f.one()));

    if (size.x() < 50) {
        shadeTriBBSingle(triData, b);
        return;
    }

    const center = bounds.center();
    // const half = bounds.halfSize();

    const dirs = [4]Vec4f{
        bounds.min,
        bounds.max,
        bounds.max.addDup(Vec4f.init(0, -size.y(), 0, 0)),
        bounds.min.addDup(Vec4f.init(0, size.y(), 0, 0)),
    };

    var totalInside: u32 = 0;

    inline for (0..dirs.len) |i| {
        const p = dirs[i];

        var triBary = Vec4f.triBarycentericCoordsOld(
            triData.screenVertex[0],
            triData.screenVertex[1],
            triData.screenVertex[2],
            p,
        );

        totalInside += @intFromBool(std.math.signbit(triBary.minElement()));
    }

    if (totalInside == dirs.len) {
        // all inside so just render the volume
        shadeTriBBSingle(triData, bounds);
        return;
    }

    inline for (0..dirs.len) |i| {
        const p = dirs[i];
        var child = Bounds.init(center, center);
        child.add(p);

        shadeTriBB(triData, child);
    }
}

const BaryEdge = struct {
    const Self = @This();
    v: Vec4f,
    nextX: Vec4f,
    nextY: Vec4f,

    pub fn init(v0: Vec4f, v1: Vec4f, start: Vec4f) Self {
        const A: f32 = v0.y() - v1.y();
        const B: f32 = v1.x() - v0.x();
        const C: f32 = v0.x() * v1.y() - v0.y() * v1.x();

        const stepX = Vec4f.splat(A * 1);
        const stepY = Vec4f.splat(B * 1);

        const x = Vec4f.splat(start.x()); // + other offsets
        const y = Vec4f.splat(start.y());

        const vA = Vec4f.splat(A);
        const vB = Vec4f.splat(B);
        const vC = Vec4f.splat(C);

        var startBary = vA.mulDup(x);
        startBary.add(vB.mulDup(y));
        startBary.add(vC);

        return .{
            .v = startBary,
            .nextX = stepX,
            .nextY = stepY,
        };
    }
};

fn shadeTriBBSingle(triData: *TriRasterData, bounds: Bounds) void {
    var srt = core.profiler.Sampler.initAndBegin(profile.?, @src().fn_name, 4);
    defer srt.end();

    var y = bounds.min.y() + 0.5;
    var p: Vec4f = Vec4f.init(0, 0, 0, 0);
    var pixelNormal: Vec4f = Vec4f.init(0, 0, 0, 0);
    var fbc: Vec4f = Vec4f.init(0, 0, 0, 1);
    var uv: Vec4f = Vec4f.init(0, 0, 0, 1);
    // const w = uv.w();
    // var c: Color = Color.black();
    const shader = triData.meshData.shader;

    _ = @atomicRmw(u32, &stats.renderedTris, .Add, 1, .seq_cst);

    const minx: i32 = @intFromFloat(bounds.min.x() + 0.5);
    const miny: i32 = @intFromFloat(bounds.min.y() + 0.5);

    var dsy = depthBuffer[currentBuffer].pxIndex(minx, miny);
    const width: usize = @intCast(colorBuffer[currentBuffer].w);
    const depthTest = shader.depthTest;

    // inline for (0..3) |i| {
    //     const pix = triData.screenVertex[i];
    //     const pixx: i32 = @intFromFloat(pix.x());
    //     const pixy: i32 = @intFromFloat(pix.y());
    //     writePixel(pixx, pixy, pix.z(), Color.fromNormalVec4f(triData.color[i]));
    // }

    const baryStart = Vec4f.triBarycentericCoordsOld(
        triData.screenVertex[0],
        triData.screenVertex[1],
        triData.screenVertex[2],
        bounds.min.addScalarDup(0.5),
    );

    // const edge12 = BaryEdge.init(triData.screenVertex[1], triData.screenVertex[2], bounds.min);
    // const edge20 = BaryEdge.init(triData.screenVertex[2], triData.screenVertex[0], bounds.min);
    // const edge01 = BaryEdge.init(triData.screenVertex[0], triData.screenVertex[1], bounds.min);

    // const edges = [3]BaryEdge{ edge12, edge20, edge01 };
    // var wrow = [3]Vec4f{ edge12.v, edge20.v, edge01.v };

    var baryRow = baryStart;
    const baryStepX = Vec4f.init(
        triData.screenVertex[2].y() - triData.screenVertex[1].y(),
        triData.screenVertex[0].y() - triData.screenVertex[2].y(),
        triData.screenVertex[1].y() - triData.screenVertex[0].y(),
        1,
    );
    const baryStepY = Vec4f.init(
        triData.screenVertex[1].x() - triData.screenVertex[2].x(),
        triData.screenVertex[2].x() - triData.screenVertex[0].x(),
        triData.screenVertex[0].x() - triData.screenVertex[1].x(),
        1,
    );
    // var minBary = Vec4f.triBarycentericCoordsOld(triData.screenVertex[0], triData.screenVertex[1], triData.screenVertex[2], bounds.min);

    while (y <= bounds.max.y()) {
        var x = bounds.min.x();
        defer y += 1;
        defer dsy += width;
        defer baryRow.add(baryStepY);

        var dsx = dsy;

        var wspan = baryRow;

        while (x <= bounds.max.x()) {
            // var sprt = core.profiler.Sampler.initAndBegin(profile.?, @src().fn_name, 2);
            // defer sprt.end();

            // const pixtracy = trace(@src());
            // defer pixtracy.end();
            // var pprof = profile.?.beginSample("render.mesh.draw.tri.pixel");
            // defer profile.?.endSample(pprof);
            //_ = @atomicRmw(u32, &stats.totalPixels, .Add, 1, .seq_cst);
            defer dsx += 1;
            defer x += 1;
            defer wspan.add(baryStepX);

            p.setX(x);
            p.setY(y);

            var triBary = wspan;

            const minBary = triBary.minElement();
            const outsideTri = std.math.signbit(minBary) or !std.math.isFinite(minBary);

            if (outsideTri) // and (triData.clipCount == 0 or (x != bounds.max.x() and x != bounds.min.x() and y != bounds.max.y() and y != bounds.min.y())))
                continue;

            triBary.div(triData.screenArea);

            const z = (triBary.x() * triData.screenVertex[0].z() +
                triBary.y() * triData.screenVertex[1].z() +
                triBary.z() * triData.screenVertex[2].z());

            p.setZ(z);

            var write = outsideTri;
            if (depthTest == 1) {
                const depth = depthBuffer[currentBuffer].readIndex(dsx);
                write = std.math.isInf(depth) or depth > z;
            }

            if (write) {

                // interpolate vertex colors across all pixels
                uv = triBary.triInterpArray(triData.uv, 0);
                const w = 1 / uv.w();

                fbc = triBary.triInterpArray(triData.color, 1);
                pixelNormal = triBary.triInterpArrayScale(triData.worldNormals, 0, w);

                const pl = triBary.triInterpArray(triData.worldVertex, 1);

                const vc = shader.pixelShader(
                    triData.meshData.model,
                    &triData.meshData.mv,
                    &triData.meshData.mvp,
                    p,
                    pl, //triData.worldNormals[0],
                    fbc,
                    pixelNormal,
                    uv.scaleDup(w),
                    shader,
                );

                triData.pixelsFilled += 1;

                const fc = Color.fromNormalVec4f(vc);
                colorBuffer[currentBuffer].writeIndex(dsx, fc);
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

        const starty = y + lines * font.glyphHeight + 2;

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
// pub fn drawPoint(mvp: *const Mat44f, point: Vec4f, color: Vec4f, shader: *Material) void {
//     const px = shader.vertexShader(mvp, 0, point, shader);
//     const pc = color;

//     const c = Color.init(@as(u8, @intFromFloat(pc.x() * 255)), @as(u8, @intFromFloat(pc.y() * 255)), @as(u8, @intFromFloat(pc.z() * 255)), @as(u8, @intFromFloat(pc.w() * 255)));

//     if (px.x() >= 0 and px.x() <= 1000 and px.y() >= 0 and px.y() <= 1000)
//         writePixel(@as(i32, @intFromFloat(px.x())), @as(i32, @intFromFloat(px.y())), c);
// }

pub fn drawPoint(pv: *const Mat44f, point: Vec4f, color: Vec4f) void {
    const spx = projectVertex(pv, point).clamped(renderBounds.min, renderBounds.max);
    // const epx = projectVertex(pv, end).clamped(renderBounds.min, renderBounds.max);

    writePixel(
        @as(i32, @intFromFloat(spx.x())),
        @as(i32, @intFromFloat(spx.y())),
        1,
        Color.fromNormalVec4f(color),
    );
}

fn projectVertex(p: *const Mat44f, v: Vec4f) Vec4f {
    const zone = trace(@src());
    defer zone.end();

    var vin = v;
    vin.setW(1);

    var out = p.mul_vec4(vin);
    out.div(out.w());

    const half = viewport.scaleDup(0.5);
    // center in viewport
    out.setX((half.x() * out.x()) + half.x());
    out.setY((half.y() * -out.y() + half.y()));
    return out;
}

///
pub fn drawWorldLine(pv: *const Mat44f, start: Vec4f, end: Vec4f, color: Vec4f) void {
    const spx = projectVertex(pv, start).clamped(renderBounds.min, renderBounds.max);
    const epx = projectVertex(pv, end).clamped(renderBounds.min, renderBounds.max);

    drawLine(
        @as(i32, @intFromFloat(spx.x())),
        @as(i32, @intFromFloat(spx.y())),
        @as(i32, @intFromFloat(epx.x())),
        @as(i32, @intFromFloat(epx.y())),
        Color.fromNormalVec4f(color),
    );
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

pub fn drawTerrain(terrain: *Terrain, model: *const Mat44f, view: *const Mat44f, proj: *const Mat44f, shader: *Material) void {

    // Model * View
    var mv = view.*;
    mv.mul(model.*);

    // View * projection
    var vp = proj.*;
    vp.mul(view.*);

    // Projection * Model * View
    var mvp = proj.*;
    mvp.mul(view.*);
    mvp.mul(model.*);

    // const size = terrain.bounds.size();
    const scale = Vec4f.init(1, 100, 1, 1);

    // const size = terrain.size;

    // const sizex: usize = @intCast(size.x());
    // const sizez: usize = 10; //@intCast(size.z());
    // for (0..sizez) |vz| {
    //     for (0..sizez) |vx| {
    //         var vscreen: [4]Vec4f = undefined;
    //         var v = [4]Vec4f{
    //             terrain.vertexBuffer[(vz + 0) * sizex + (vx + 0)],
    //             terrain.vertexBuffer[(vz + 1) * sizex + (vx + 0)],
    //             terrain.vertexBuffer[(vz + 1) * sizex + (vx + 1)],
    //             terrain.vertexBuffer[(vz + 0) * sizex + (vx + 1)],
    //         };
    //         var outsideCount: usize = 0;
    //         var screenBounds = Bounds.initInfinity();
    //         inline for (0..v.len) |i| {
    //             const cameraVertex = shader.vertexShader(&mv, @intCast(i), v[i].mulDup(scale), shader);
    //             const clipVertex = shader.projectionShader(proj, cameraVertex, viewport, shader);
    //             var screenVertex = clipVertex.divDup(clipVertex.w());

    //             if (screenVertex.minElement() < -1 or screenVertex.maxElement() > 1)
    //                 outsideCount += 1;

    //             const half = viewport.scaleDup(0.5);
    //             const out = screenVertex;
    //             screenVertex.setX(half.x() * out.x() + half.x());
    //             screenVertex.setY(half.y() * -out.y() + half.y());

    //             vscreen[i] = screenVertex;
    //             screenBounds.add(screenVertex);
    //         }

    //         // if (outsideCount == 4)
    //         //     continue;

    //         screenBounds.limit(renderBounds);
    //         screenBounds.topLeftHandLimit();

    //         for (0..10) |i| {
    //             const ioffset: f32 = @as(f32, @floatFromInt(i)) / 10.0;
    //             for (0..10) |j| {
    //                 const joffset = @as(f32, @floatFromInt(j)) / 10.0;
    //                 const offset = Vec4f.init(ioffset, 0, joffset, 0);
    //                 const p = v[0].addDup(offset);

    //                 var h: f32 = 0;
    //                 var hsum: f32 = 0;
    //                 inline for (0..v.len) |vi| {
    //                     var dir = v[vi].subDup(p);
    //                     dir.setY(0);
    //                     const dr = dir.length3();
    //                     const d = @min(1.0, 1.0 / dr);
    //                     h += d * v[vi].y();
    //                     hsum += v[vi].y();
    //                 }

    //                 h /= hsum;
    //                 const pfinal = Vec4f.init(p.x(), h * scale.y(), p.z(), 1);

    //                 const cameraVertex = shader.vertexShader(&mv, 0, pfinal, shader);
    //                 const clipVertex = shader.projectionShader(proj, cameraVertex, viewport, shader);
    //                 var screenVertex = clipVertex.divDup(clipVertex.w());

    //                 if (screenVertex.minElement() < -1 or screenVertex.maxElement() > 1)
    //                     continue;

    //                 const half = viewport.scaleDup(0.5);
    //                 const out = screenVertex;
    //                 screenVertex.setX(half.x() * out.x() + half.x());
    //                 screenVertex.setY(half.y() * -out.y() + half.y());

    //                 const fc = Color.fromNormalVec4f(Vec4f.splat3((h) * 2, 1));
    //                 const pixel = core.vector.Vec4fToVec4i(screenVertex);

    //                 colorBuffer[currentBuffer].write(pixel.x(), pixel.y(), fc);
    //             }
    //         }
    //     }
    // }

    for (0..terrain.vertexBuffer.len) |i| {
        const v = terrain.vertexBuffer[i];
        //const vn = terrain.vertexBuffer[i + 1];
        const cameraVertex = shader.vertexShader(&mv, @intCast(i), v.mulDup(scale), shader);

        const clipVertex = shader.projectionShader(proj, cameraVertex, viewport, shader);

        var screenVertex = clipVertex.divDup(clipVertex.w());

        if (screenVertex.minElement() < -1 or screenVertex.maxElement() > 1)
            continue;

        const half = viewport.scaleDup(0.5);
        const out = screenVertex;
        screenVertex.setX(half.x() * out.x() + half.x());
        screenVertex.setY(half.y() * -out.y() + half.y());

        // const c = shader.pixelShader(
        //     model,
        //     &mv,
        //     &mvp,
        //     screenVertex,
        //     screenVertex,
        //     v,
        //     Vec4f.up(),
        //     Vec4f.init(v.x() / size.x(), v.z() / size.z(), 0, 1),
        //     shader,
        // );

        const fc = Color.fromNormalVec4f(Vec4f.splat3(v.y() * 2, 1));
        const pixel = core.vector.Vec4fToVec4i(screenVertex);

        colorBuffer[currentBuffer].write(pixel.x(), pixel.y(), fc);
        // depthBuffer[currentBuffer].writeIndex(dsx, z);
    }
}
