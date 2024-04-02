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

const matfuncs = engine.render.material;

var modelMat = engine.Mat44f.identity();
var viewMat = engine.Mat44f.identity();
var mesh: engine.Mesh = undefined;
var projMat: engine.Mat44f = undefined;

var meshAllocator = std.heap.page_allocator;
var textureAllocator = std.heap.page_allocator;
var meshMaterial: engine.render.Material = undefined;

var render3d: bool = true;
var renderSingleFrame: bool = false;

var nearZ: f32 = 1;
var farZ: f32 = 10000.0;
var fovY: f32 = 90.0;
var aspect: f32 = 0.0;

pub fn init() !void {
    aspect = @as(f32, @floatFromInt(engine.systemConfig.renderWidth)) / @as(f32, @floatFromInt(engine.systemConfig.renderHeight));
    projMat = engine.Mat44f.createPerspective(fovY, aspect, nearZ, farZ);
    // projMat = engine.Mat44f.createPerspectiveSimple(fovY, aspect, nearZ, farZ);

    // mesh = try tools.MeshObjLoader.importObjFile(&meshAllocator, "../../assets/cube.obj");
    // mesh = try tools.MeshObjLoader.importObjFile(&meshAllocator, "../../assets/triangle.obj");
    // mesh = try tools.MeshObjLoader.importObjFile(&meshAllocator, "../../assets/bed.obj");
    //mesh = try tools.MeshObjLoader.importObjFile(&meshAllocator, "../../assets/crates/crate-04-1.obj");
    // mesh = try tools.MeshObjLoader.importObjFile(&meshAllocator, "../../assets/suzanne.obj");
    mesh = try tools.MeshObjLoader.importObjFile(&meshAllocator, "../../assets/E1M1.bsp.geometry.tri.obj");
    //mesh = try tools.MeshObjLoader.importObjFile(&meshAllocator, "../../assets/plane.obj");
    // mesh = try tools.MeshObjLoader.importObjFile(&meshAllocator, "../../assets/axis.obj");

    //mesh = try tools.MeshObjLoader.importObjFile(&meshAllocator, "../../assets/Character.obj");
    // var  texture = try tools.TgaTexLoader.importTGAFile(textureAllocator, "../../assets/black_rock.tga");
    var texture = try tools.TgaTexLoader.importTGAFile(&textureAllocator, "../../assets/grass.tga");
    // var texture = try tools.TgaTexLoader.importTGAFile(&textureAllocator, "../../assets/crates/diffuse.tga");

    meshMaterial = engine.render.Material{
        .backfaceCull = true,
        .depthTest = 1,
        .lightDirection = engine.Vec4f.init(-0.913913, 0.389759, -0.113369, 1).normalized3(),
        .lightColor = engine.Vec4f.one(),
        .lightIntensity = 1.2,
        .vertexShader = applyVertexShader,
        .projectionShader = projectVertexNoDiv,
        .pixelShader = shader_lit_texture0_frag,
        .texture = texture,
    };

    var fontTex = try tools.TgaTexLoader.importTGAFile(&textureAllocator, "../../assets/mbf_small_7x7.tga");

    font = engine.render.Font{
        .glyphWidth = 7,
        .glyphHeight = 7,
        .texture = fontTex,
    };

    viewMat.translate(engine.Vec4f.init(0, 0, -10.0, 0));
}

pub fn shutdown() !void {
    textureAllocator.deinit();
    meshAllocator.deinit();
}

const moveSpeed = 10;

var mousePos = engine.Vec4f.zero();
var cameraPos = engine.Vec4f.zero();
var cameraRot = engine.Mat44f.identity();

var exposure_bias: f32 = 2.0;
var font: engine.render.Font = undefined;
var singleFrameKeyDown: bool = false;

pub fn update() bool {
    // const tracy = trace(@src());
    // defer tracy.end();

    if (input.isKeyDown(input.KeyCode.ESCAPE))
        return false;

    var currentMouse = engine.Vec4f.init(((@as(f32, @floatFromInt(input.getMouseX())) / @as(f32, @floatFromInt(engine.systemConfig.windowWidth))) - 0.5) * 2.0, ((@as(f32, @floatFromInt(input.getMouseY())) / @as(f32, @floatFromInt(engine.systemConfig.windowHeight))) - 0.5) * 2.0, 0, 0);

    const mouseDelta = currentMouse.subDup(mousePos);

    const depth = (input.keyStateFloat(input.KeyCode.W) - input.keyStateFloat(input.KeyCode.S)) * moveSpeed;
    const horizontal = (input.keyStateFloat(input.KeyCode.D) - input.keyStateFloat(input.KeyCode.A)) * moveSpeed;
    const vertical = (input.keyStateFloat(input.KeyCode.UP) - input.keyStateFloat(input.KeyCode.DOWN)) * moveSpeed;

    // const rot = (input.keyStateFloat(input.KeyCode.Q) - input.keyStateFloat(input.KeyCode.E)) * moveSpeed;

    const rightButton = @as(f32, @floatFromInt(input.getMouseRight()));

    const yaw = mouseDelta.x();
    const pitch = mouseDelta.y();

    currentMouse.setZ(yaw);
    currentMouse.setW(pitch);

    const exposure = (input.keyStateFloat(input.KeyCode.U) - input.keyStateFloat(input.KeyCode.J)) * 0.1;
    const bright = input.keyStateFloat(input.KeyCode.I) - input.keyStateFloat(input.KeyCode.K) * 0.1;
    meshMaterial.lightIntensity = @max(meshMaterial.lightIntensity + bright, 0.0);
    exposure_bias = @max(exposure_bias + exposure, 0.0);

    mousePos = currentMouse;
    var trans = engine.Mat44f.identity();

    //trans.mul33(viewMat);
    trans.translate(engine.Vec4f.init(-horizontal, -vertical, depth, 0));
    trans.mul(engine.Mat44f.rotX(rightButton * pitch));
    trans.mul(engine.Mat44f.rotY(rightButton * yaw));
    trans.mul(viewMat);

    viewMat = trans;

    //trans.translate(viewMat.forward().scaleDup(10));
    engine.render.viewFrustum.from(
        &trans,
        aspect,
        fovY,
        nearZ,
        farZ,
    );

    // var rotmat = engine.Mat44f.rotY(0.1 / 60.0);

    // modelMat.mul(rotmat);

    _ = engine.sys.showMouseCursor(~input.getMouseRight());
    _ = engine.sys.setRelativeMouseMode(input.getMouseRight());

    // var srenderDraw = engine.Sampler.begin(&engine.profiler,"draw.mesh");
    // defer srenderDraw.end();

    if (!input.isKeyDown(input.KeyCode.SPACE)) {
        // if(renderSingleFrame)
        //     render3d = false;

        // const renderStart = frameTimer.read();
        // renderTimer.reset();
        engine.render.drawMesh(&modelMat, &viewMat, &projMat, &mesh, &meshMaterial) catch {};
        // engine.render.drawMesh(&modelMat, &viewMat, &projMat, &mesh, &meshMaterial);
        // engine.render.drawMesh(&modelMat, &viewMat, &projMat, &mesh, &meshMaterial);
        // engine.render.drawMesh(&modelMat, &viewMat, &projMat, &mesh, &meshMaterial);
        // engine.render.drawMesh(&modelMat, &viewMat, &projMat, &mesh, &meshMaterial);
        // engine.render.drawMesh(&modelMat, &viewMat, &projMat, &mesh, &meshMaterial);
        // engine.render.drawMesh(&modelMat, &viewMat, &projMat, &mesh, &meshMaterial);

        engine.render.drawString(&font, "Hello World!", 10, 10, engine.Vec4f.one());
    }

    return true;
}

fn projectVertex(p: *const engine.Mat44f, v: engine.Vec4f, viewport: engine.Vec4f, material: *engine.render.Material) engine.Vec4f {
    const zone = trace(@src());
    defer zone.end();

    var out = p.mul33_divW_vec4(v);

    const half = viewport.scaleDup(0.5);
    _ = material;
    // _ = viewport;
    // center in viewport
    out.setX((half.x() * out.x()) + half.x());
    out.setY((half.y() * -out.y() + half.y()));
    return out;
}

fn projectVertexNoDiv(
    p: *const engine.Mat44f,
    v: engine.Vec4f,
    viewport: engine.Vec4f,
    material: *engine.render.Material,
) engine.Vec4f {
    var out = p.mul_vec4(v);
    // if (v.w() != 0.0)
    //     out.div(v.w());

    // const half = viewport.scaleDup(0.5);
    _ = material;
    _ = viewport;
    // // // center in viewport
    // out.setX(half.x() * out.x() + half.x());
    // out.setY(half.y() * -out.y() + half.y());
    return out;
}

fn applyVertexShader(mv: *const engine.Mat44f, index: u16, vertex: engine.Vec4f, material: *engine.render.Material) engine.Vec4f {
    const zone = trace(@src());
    defer zone.end();

    const out = mv.mul_vec4(vertex);

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

    const W = engine.Vec4f.init(11.2, 11.2, 11.2, 0);
    const white_scale = engine.Vec4f.one().divVecDup(uncharted2_tonemap_partial(W));
    return curr.mulDup(white_scale);
}

pub inline fn reinhard(c: engine.Vec4f) engine.Vec4f {
    return c.divVecDup(engine.Vec4f.one().addDup(c));
    //return v / (1.0f + v);
}

///
fn shader_unlit_colors_frag(
    mv: *const engine.Mat44f,
    mvp: *const engine.Mat44f,
    pixel: engine.Vec4f,
    color: engine.Vec4f,
    normal: engine.Vec4f,
    uv: engine.Vec4f,
    material: *engine.render.Material,
) engine.Vec4f {
    _ = pixel;
    _ = mvp;
    _ = uv;
    _ = mv;
    _ = material;
    _ = normal;

    var c = color;
    return c;
}

fn shader_lit_texture0_frag(
    mv: *const engine.Mat44f,
    mvp: *const engine.Mat44f,
    pixel: engine.Vec4f,
    color: engine.Vec4f,
    normal: engine.Vec4f,
    uv: engine.Vec4f,
    material: *engine.render.Material,
) engine.Vec4f {

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
    // var c = material.texture.sample(uv.x(), uv.y());
    // var c = color;

    c.lerp(color, 0.25);
    c.setW(1.0);

    const l = @max(normal.dot3(material.lightDirection) * material.lightIntensity, 0.5);
    c.scale(l);
    c.setW(1.0);

    // return c;
    // return uncharted2_filmic(c);
    return reinhard(c);
    //return c.scaleDup(l);
}
