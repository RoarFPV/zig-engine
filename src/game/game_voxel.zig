// platform imports
const std = @import("std");
const fmt = std.fmt;
const warn = std.debug.print;
const assert = std.debug.assert;

// engine imports
const engine = @import("../engine.zig");
const tools = @import("../tools.zig");
const input = engine.input;

pub const trace = @import("../tracy.zig").trace;
pub const World = @import("world.zig").World;
const matfuncs = engine.render.material;

const sdf = @import("../core/sdf.zig");

var modelMat = engine.Mat44f.identity();
var viewMat = engine.Mat44f.identity();
var mesh: engine.Mesh = undefined;
var projMat: engine.Mat44f = undefined;

const Vec4f = engine.Vec4f;
const Mat44f = engine.Mat44f;
const Color = engine.Color;
const math = std.math;

var worldAllocator = std.heap.page_allocator;
var textureAllocator = std.heap.page_allocator;
var meshMaterial: engine.render.Material = undefined;

var render3d: bool = true;
var renderSingleFrame: bool = false;

var world: World = undefined;
var mainTimer: std.time.Timer = undefined;

pub fn init() !void {
    projMat = engine.Mat44f.createPerspective(50, @as(f32, @floatFromInt(engine.systemConfig.renderWidth)) / @as(f32, @floatFromInt(engine.systemConfig.renderHeight)), 0.1, 1000);

    meshMaterial = engine.render.Material{
        .depthTest = 1,
        .lightDirection = engine.Vec4f.init(-0.913913, 0.389759, -0.113369, 1).normalized3(),
        .lightColor = engine.Vec4f.one(),
        .lightIntensity = 1,
        .vertexShader = applyVertexShader,
        .projectionShader = projectVertex,
        .pixelShader = applyPixelShader,
        .texture = undefined,
        .viewPos = Vec4f.zero(),
        .roughness = 0.5,
        .metal = 0.2,
    };

    var fontTex = try tools.TgaTexLoader.importTGAFile(&textureAllocator, "../../assets/mbf_small_7x7.tga");

    font = engine.render.Font{
        .glyphWidth = 7,
        .glyphHeight = 7,
        .texture = fontTex,
    };

    viewMat.translate(engine.Vec4f.init(0, 1, -1.0, 0));

    world = World.init(worldAllocator);

    mainTimer = try std.time.Timer.start();

    const size = 2;
    for (0..size) |y| {
        for (0..size) |x| {
            for (0..size) |z| {
                try world.nodes.put(@Vector(4, i32){ @as(i32, @intCast(x)), @as(i32, @intCast(y)), @as(i32, @intCast(z)), 1 }, @as(u32, @intCast(y * size + x)));
            }
        }
    }
}

pub fn shutdown() !void {
    world.deinit();
    textureAllocator.deinit();
    worldAllocator.deinit();
}

const moveSpeed = 0.1;

var mousePos = engine.Vec4f.zero();
var cameraPos = engine.Vec4f.zero();
var cameraRot = engine.Mat44f.identity();

var exposure_bias: f32 = 2.0;
var font: engine.render.Font = undefined;
var singleFrameKeyDown: bool = false;

var runTime: f32 = 0.0;

pub fn update() bool {
    // const tracy = trace(@src());
    // defer tracy.end();

    runTime = @as(f32, @floatCast(@as(f64, @floatFromInt(mainTimer.read())) / 1_000_000_000.0));

    if (input.isKeyDown(input.KeyCode.ESCAPE))
        return false;

    var currentMouse = engine.Vec4f.init(((@as(f32, @floatFromInt(input.getMouseX())) / @as(f32, @floatFromInt(engine.systemConfig.windowWidth))) - 0.5) * 2.0, ((@as(f32, @floatFromInt(input.getMouseY())) / @as(f32, @floatFromInt(engine.systemConfig.windowHeight))) - 0.5) * 2.0, 0, 0);

    const mouseDelta = currentMouse.subDup(mousePos);
    const rightButton = @as(f32, @floatFromInt(input.getMouseRight()));

    var forward: f32 = moveSpeed;
    if (rightButton > 0.0)
        forward *= 0.1;

    const depth = (input.keyStateFloat(input.KeyCode.W) - input.keyStateFloat(input.KeyCode.S)) * forward;
    const horizontal = (input.keyStateFloat(input.KeyCode.D) - input.keyStateFloat(input.KeyCode.A)) * moveSpeed;
    const vertical = (input.keyStateFloat(input.KeyCode.UP) - input.keyStateFloat(input.KeyCode.DOWN)) * moveSpeed;

    // const rot = (input.keyStateFloat(input.KeyCode.Q) - input.keyStateFloat(input.KeyCode.E)) * moveSpeed;

    const yaw = mouseDelta.x();
    const pitch = mouseDelta.y();

    currentMouse.setZ(yaw);
    currentMouse.setW(pitch);

    const exposure = (input.keyStateFloat(input.KeyCode.U) - input.keyStateFloat(input.KeyCode.J)) * moveSpeed;
    const bright = input.keyStateFloat(input.KeyCode.I) - input.keyStateFloat(input.KeyCode.K) * 0.1;
    meshMaterial.lightIntensity = @max(meshMaterial.lightIntensity + bright, 0.0);
    exposure_bias = @max(exposure_bias + exposure, 0.0);

    mousePos = currentMouse;
    var trans = engine.Mat44f.identity();

    //trans.mul33(viewMat);
    trans.translate(engine.Vec4f.init(horizontal, vertical, depth, 0));
    trans.mul(engine.Mat44f.rotX(rightButton * pitch));
    trans.mul(engine.Mat44f.rotY(rightButton * yaw));
    trans.mul(viewMat);

    viewMat = trans;

    // var pos = modelMat.
    //var rotmat = engine.Mat44f.rotY(0.01 / 60.0);

    //modelMat.mul(rotmat);

    //modelMat = rotmat;

    _ = engine.sys.showMouseCursor(~input.getMouseRight());
    _ = engine.sys.setRelativeMouseMode(input.getMouseRight());

    // var srenderDraw = engine.Sampler.begin(&engine.profiler,"draw.mesh");
    // defer srenderDraw.end();

    meshMaterial.lightDirection = Vec4f.init(lightPos.x() - 2 * @sin(runTime + 2), 10.0, lightPos.x() + 2 * @cos(runTime), 0);
    lightPos = meshMaterial.lightDirection;
    var buffer: [100]u8 = undefined;
    const buf = buffer[0..];

    if (!input.isKeyDown(input.KeyCode.SPACE)) {
        //renderWorld(projMat, viewMat.mulDup(modelMat));

        engine.render.renderFrame(&viewMat, rayMarchPixel);

        const pos = viewMat.position();
        const s = std.fmt.bufPrint(buf, "o: {0:.3}\n{1:.3}\n{2:.3}", .{ pos.x(), pos.y(), pos.z() }) catch unreachable;

        engine.render.drawString(&font, s, 10, 10, engine.Vec4f.one());
    }

    return true;
}

pub fn map(p: Vec4f) f32 {
    //const pivot = p.addDup(Vec4f.init(0,0,10, 0));
    //var d = p.y();//sdf.Plane(p, Vec4f.init(0,1,0,0).normalized(), 0);
    var d: f32 = maxDist + 1;
    d = @min(d, sdf.Box(p.subDup(Vec4f.init(0, -2, 0, 0)), Vec4f.init(100, 0.1, 100, 0)));
    d = @min(d, sdf.Cube(p.subDup(Vec4f.init(-2, 1, 6, 0)), 0.25));
    return @min(d, sdf.Sphere(p.subDup(Vec4f.init(0, 1, 6, 0)), 1));
}

pub fn march(origin: Vec4f, ray: Vec4f, maxSteps: usize) f32 {
    var dist: f32 = 0.0;

    for (0..maxSteps) |s| {
        _ = s;
        const p = origin.addDup(ray.scaleDup(dist));
        const d = map(p);
        dist += d;

        if (d < minSurfaceDist)
            return dist;

        if (dist >= maxDist)
            return maxDist;
    }

    return dist;
}

var lightPos = Vec4f.init(10, -1, 4, 0);
const MaxSteps: usize = 50;
const maxDist: f32 = 1000.0;
const minSurfaceDist = 0.001;

pub fn lightAt(p: Vec4f, n: Vec4f, light: Vec4f) f32 {
    _ = light;

    const material = meshMaterial;
    //const lr = light.subDup(p);
    //const l = lightPos.normalized();// lr.normalized();

    //const diff = math.clamp(n.dot(l), 0.1, 1.0);

    var diff = @max(n.dot3(material.lightDirection) * material.lightIntensity, 0.4);
    //c.scale(diff);

    const d = march(p.addDup(n.scaleDup(minSurfaceDist * 3)), material.lightDirection, 10);

    if (d < maxDist)
        diff *= 0.1;

    return diff;
}

pub fn rayMarchPixel(view: *const Mat44f, fragCoord: Vec4f, viewport: Vec4f) Color {
    //_=view;
    //var rayOrigin = Vec4f.init(0,1,0,0);
    var rayOrigin = view.position();

    // // center in viewport
    const ray = Vec4f.init((fragCoord.x() - viewport.x() * 0.5) / viewport.y(), (viewport.y() * 0.5 - fragCoord.y()) / viewport.y(), 1, 0).normalized();

    const dist = march(rayOrigin, ray, MaxSteps);

    if (dist >= maxDist)
        return Color.black();

    const p = rayOrigin.addDup(ray.scaleDup(dist));
    const n = sdf.normal2(p, dist, map);

    const c = lightAt(p, n, lightPos);

    const vc = Vec4f.splat(c);

    return Color.fromNormalVec4f(uncharted2_filmic(vc));
}

pub fn renderWorld(proj: engine.Mat44f, mv: engine.Mat44f) void {
    var it = world.nodes.iterator();
    var count: usize = 0;
    while (it.next()) |entry| : (count += 1) {
        const key = entry.key_ptr.*;

        renderNode(&proj, &mv, key, entry.value_ptr.*);
    }
}

pub fn renderNode(proj: *const engine.Mat44f, mv: *const engine.Mat44f, loc: @Vector(4, i32), node: u32) void {
    if (node <= 0)
        return;

    const p = engine.Vec4f.init(@as(f32, @floatFromInt(loc[0])), @as(f32, @floatFromInt(loc[1])), @as(f32, @floatFromInt(loc[2])), @as(f32, @floatFromInt(loc[3])));

    engine.render.drawSdfBox(proj, mv, p, engine.Vec4f.init(p.z() / 10, p.y() / 100, p.x() / 10, 1.0), &meshMaterial, Vec4f.one());
}

fn projectVertex(p: *const engine.Mat44f, v: engine.Vec4f, viewport: engine.Vec4f, material: *engine.render.Material) engine.Vec4f {
    var out = p.mul33_vec4(v);
    const half = viewport.scaleDup(0.5);
    _ = material;

    // center in viewport
    out.setX(half.x() * out.x() + half.x());
    out.setY(half.y() * -out.y() + half.y());
    return out;
}

fn applyVertexShader(mvp: *const engine.Mat44f, index: u16, vertex: engine.Vec4f, material: *engine.render.Material) engine.Vec4f {
    const out = mvp.mul_vec4(vertex);
    _ = material;
    _ = index;
    return out;
}

pub inline fn uncharted2_tonemap_partial(x: engine.Vec4f) engine.Vec4f {
    const A = 0.15;
    const B = 0.50;
    const C = 0.10;
    const D = 0.20;
    const E = 0.02;
    const F = 0.30;

    const EdivF = E / F;
    const DmulE = D * E;
    const DmulF = D * F;
    const CmulB = C * B;

    const xmulA = x.scaleDup(A);

    var xNumer = x.mulDup(xmulA.addScalarDup(CmulB));
    xNumer.addScalar(DmulE);

    var xDenom = x.mulDup(xmulA.addScalarDup(B));
    xDenom.addScalar(DmulF);

    xNumer.divVec(xDenom);
    xNumer.subScalar(EdivF);

    return xNumer;
    //return ((x*(A*x+C*B)+D*E)/(x*(A*x+B)+D*F))-E/F;
}

pub inline fn uncharted2_filmic(v: engine.Vec4f) engine.Vec4f {
    //const exposure_bias = 2.0;
    const curr = uncharted2_tonemap_partial(v.scaleDup(exposure_bias));

    const W = engine.Vec4f.init(11.2, 11.2, 11.2, 1);
    const tmp = uncharted2_tonemap_partial(W);
    const white_scale = engine.Vec4f.one().divVecDup(tmp);
    return curr.mulDup(white_scale);
}

pub inline fn reinhard(c: engine.Vec4f) engine.Vec4f {
    return c.divVecDup(engine.Vec4f.one().addDup(c));
    //return v / (1.0f + v);
}

///
fn applyPixelShader(
    m: *const Mat44f,
    mv: *const Mat44f,
    mvp: *const Mat44f,
    pixel: Vec4f,
    pixelLocal: Vec4f,
    color: Vec4f,
    normal: Vec4f,
    uv: Vec4f,
    material: *engine.render.Material,
) engine.Vec4f {
    // var c = color.addDup(
    //   engine.Vec4f.init(
    //     (std.math.sin(uv.x*uv.y*1000)+1/2),
    //     (std.math.cos(uv.y*1000)+1/2),
    //     0,1)
    //   );

    _ = m;
    _ = pixelLocal;
    _ = mv;
    _ = color;
    _ = pixel;
    _ = mvp;
    //var c = material.texture.sample(uv.x(), uv.y());
    var c = material.texture.sampleBilinear(uv.x(), uv.y());

    const l = @max(normal.dot3(material.lightDirection) * material.lightIntensity, 0.4);

    c.scale(l);

    //return c;
    return uncharted2_filmic(c);
    //kreturn reinhard(c);
    //return c.scaleDup(l);
}
