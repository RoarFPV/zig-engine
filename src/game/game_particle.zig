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
var particleList: ParticleList = undefined;
const cubeSize = 5;
var particleCount: usize = cubeSize * cubeSize * cubeSize;
var particleWorld = ParticleWorld{};

const smoothingRadius: f32 = 4;
const smoothingRadiusSqr: f32 = smoothingRadius * smoothingRadius;

const worldSize = 6;
const worldSizeHalf = worldSize / 2;
const boundDamping = -0.8;

const ParticleWorld = struct {
    const Bounds = struct {
        min: Vec4f = Vec4f.init(-worldSizeHalf, -worldSizeHalf, -worldSizeHalf / 2, 1),
        max: Vec4f = Vec4f.init(worldSizeHalf, worldSizeHalf, worldSizeHalf / 2, 1),
    };

    // smoothing kernels defined in Müller and their gradients
    // adapted to 2D per "SPH Based Shallow Water Simulation" by Solenthaler et al.
    POLY6: f32 = 4.0 / (std.math.pi * std.math.pow(f32, smoothingRadius, 8.0)),
    SPIKY_GRAD: f32 = -10.0 / (std.math.pi * std.math.pow(f32, smoothingRadius, 5.0)),
    VISC_LAP: f32 = 4000.0 / (std.math.pi * std.math.pow(f32, smoothingRadius, 5.0)),
    VISC: f32 = 0.005,

    bounds: Bounds = Bounds{},
    gravity: Vec4f = Vec4f.init(0, -10.0, 0, 0),
    gas: f32 = 10.0,
    restingDensity: f32 = 0, //30.0, //300,
};

const Particle = struct {
    index: usize,
    position: Vec4f = Vec4f.zero(),
    predictedPos: Vec4f = Vec4f.zero(),
    gridPos: Vec4f = Vec4f.zero(),

    velocity: Vec4f = Vec4f.zero(),
    acceleration: Vec4f = Vec4f.zero(),
    color: Vec4f,

    density: f32 = 0.0,
    pressure: f32 = 0.0,
    mass: f32 = 1.0,
    size: f32 = 1.0,

    pub fn init(index: usize) Particle {
        return .{
            .index = index,
            .color = Vec4f.init(1, 1, 1, 1),
            .position = Vec4f.init(@floatFromInt(index), 1, 0, 1),
            .acceleration = Vec4f.zero(),
            .velocity = Vec4f.zero(),
            .gridPos = Vec4f.zero(),
        };
    }
};

const ParticleList = std.MultiArrayList(Particle);

const moveSpeed = 5;

var mousePos = engine.Vec4f.zero();
var cameraPos = engine.Vec4f.zero();
var cameraRot = engine.Mat44f.identity();

var exposure_bias: f32 = 2.0;
var font: engine.render.Font = undefined;
var singleFrameKeyDown: bool = false;

var runTime: f32 = 0.0;
var simulationStep: bool = false;

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
        .roughness = 0.4,
        .metal = 0.1,
    };

    var fontTex = try tools.TgaTexLoader.importTGAFile(&textureAllocator, "../../assets/mbf_small_7x7.tga");

    font = engine.render.Font{
        .glyphWidth = 7,
        .glyphHeight = 7,
        .texture = fontTex,
    };

    viewMat.translate(engine.Vec4f.init(0, 1, -20.0, 0));

    world = World.init(worldAllocator);

    mainTimer = try std.time.Timer.start();

    particleList = ParticleList{};
    try particleList.ensureTotalCapacity(worldAllocator, particleCount);

    var i: u32 = 0;
    for (0..cubeSize) |z| {
        for (0..cubeSize) |y| {
            for (0..cubeSize) |x| {
                var p = Particle.init(i);
                p.position = Vec4f.init(
                    @as(f32, @floatFromInt(x)) * 0.5,
                    @as(f32, @floatFromInt(y)) * 0.5,
                    @as(f32, @floatFromInt(z)) * 0.5,
                    1.0,
                );
                particleList.appendAssumeCapacity(p);
                i += 1;
            }
        }
    }

    resetParticles(particleList);
    // const size = 2;
    // for (0..size) |y| {
    //     for (0..size) |x| {
    //         for (0..size) |z| {
    //             try world.nodes.put(@Vector(4, i32){ @as(i32, @intCast(x)), @as(i32, @intCast(y)), @as(i32, @intCast(z)), 1 }, @as(u32, @intCast(y * size + x)));
    //         }
    //     }
    // }
}

pub fn shutdown() !void {
    world.deinit();
    textureAllocator.deinit();
    worldAllocator.deinit();
}

pub fn update() bool {
    // const tracy = trace(@src());
    // defer tracy.end();
    var newTime = @as(f32, @floatCast(@as(f64, @floatFromInt(mainTimer.read())) / 1_000_000_000.0));
    var dt = newTime - runTime;
    runTime = newTime;

    if (input.isKeyDown(input.KeyCode.ESCAPE))
        return false;

    var currentMouse = engine.Vec4f.init(((@as(f32, @floatFromInt(input.getMouseX())) / @as(f32, @floatFromInt(engine.systemConfig.windowWidth))) - 0.5) * 2.0, ((@as(f32, @floatFromInt(input.getMouseY())) / @as(f32, @floatFromInt(engine.systemConfig.windowHeight))) - 0.5) * 2.0, 0, 0);

    const mouseDelta = currentMouse.subDup(mousePos);
    const rightButton = @as(f32, @floatFromInt(input.getMouseRight()));

    var forward: f32 = moveSpeed;
    // if (rightButton > 0.0)
    //     forward *= 0.1;

    const depth = (input.keyStateFloat(input.KeyCode.W) - input.keyStateFloat(input.KeyCode.S)) * forward;
    const horizontal = (input.keyStateFloat(input.KeyCode.A) - input.keyStateFloat(input.KeyCode.D)) * moveSpeed;
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
    trans.translate(engine.Vec4f.init(horizontal, vertical, depth, 1).scale3Dup(dt));
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

        engine.render.renderFrame(&viewMat, renderPixel);

        if (input.isKeyDown(input.KeyCode.R)) {
            resetParticles(particleList);
        }

        const gaspress = (input.keyStateFloat(input.KeyCode.P) - input.keyStateFloat(input.KeyCode.SEMICOLON));
        particleWorld.gas += gaspress * dt * 10;

        const visc = (input.keyStateFloat(input.KeyCode.O) - input.keyStateFloat(input.KeyCode.L));

        particleWorld.VISC += visc * dt * 0.1;

        const stepTime = 1.0 / 100.0;

        var remaining = @min(dt, 1.0 / 60.0);

        while (remaining > 0) {

            // if (input.isKeyDown(input.KeyCode.F)) {
            //     if (!simulationStep) {
            predict(particleList, stepTime);
            updateDensity(particleList, particleWorld);
            updatePressureForces(particleList, particleWorld);
            simulateParticles(particleList, stepTime);

            if (remaining > stepTime) {
                remaining -= stepTime;
            } else {
                remaining = 0;
            }
        }
        //         simulationStep = true;
        //     }
        // } else {
        //     simulationStep = false;
        // }

        renderParticles(particleList, dt);

        const pos = viewMat.position();
        const s = std.fmt.bufPrint(buf, "o: {0:.3}\n{1:.3}\n{2:.3}", .{ pos.x(), pos.y(), pos.z() }) catch unreachable;

        engine.render.drawString(&font, s, 10, 10, engine.Vec4f.one());
    }

    return true;
}

fn renderParticles(particles: ParticleList, dt: f32) void {
    _ = dt;
    var position = particles.items(.position);
    // var color = particles.items(.color);
    var size = particles.items(.size);
    var velocity = particles.items(.velocity);
    var acceleration = particles.items(.acceleration);

    for (0..particles.len, size, position, velocity, acceleration) |p, s, pos, vel, acc| {
        _ = acc;
        _ = p;
        _ = s;

        const speedSq = vel.length3Sqr() / 100.0;
        const color = Vec4f.init(
            @mod(speedSq, 0.75),
            @mod(speedSq, 0.5),
            0.5 + @mod(speedSq, 0.25),
            1.0,
        );

        // engine.render
        //     .drawWorldLine(&projMat, &viewMat, pos, pos.addDup(vel)
        //     .scaleDup(1 / vel.length3Sqr()), Vec4f.init(1.0, 0, 0, 1.0), &meshMaterial);

        // engine.render
        //     .drawWorldLine(&projMat, &viewMat, pos, pos.addDup(acc)
        //     .scaleDup(1), Vec4f.init(0.0, 1.0, 0, 1.0), &meshMaterial);

        engine.render.drawPoint(&projMat, &viewMat, pos, color, &meshMaterial);
    }
}

fn resetParticles(particles: ParticleList) void {
    var i: u32 = 0;
    var pos = particles.items(.position);
    var acc = particles.items(.acceleration);
    var vel = particles.items(.velocity);
    var dens = particles.items(.density);
    var press = particles.items(.pressure);

    for (0..cubeSize) |z| {
        for (0..cubeSize) |y| {
            for (0..cubeSize) |x| {
                var wy: f32 = @floatFromInt(y);
                var wx: f32 = @floatFromInt(x);
                var wz: f32 = @floatFromInt(z);
                pos[i] = Vec4f.init(wx * 0.5, wy * 0.5, wz * 0.5, 1.0);
                vel[i] = Vec4f.zero();
                acc[i] = Vec4f.zero();
                dens[i] = 0.5;
                press[i] = 0.0;
                i += 1;
            }
        }
    }
}

fn simulateParticles(particles: ParticleList, dt: f32) void {
    var position = particles.items(.position);
    var predicted = particles.items(.predictedPos);
    var velocity = particles.items(.velocity);
    var acceleration = particles.items(.acceleration);
    var mass = particles.items(.mass);
    var density = particles.items(.density);
    const bounds = particleWorld.bounds;

    for (
        0..particles.len,
        acceleration,
        velocity,
        position,
        predicted,
        mass,
        density,
    ) |p, *acc, *vel, *pos, *pred, m, dens| {
        _ = pred;
        _ = p;

        var grav = particleWorld.gravity;
        grav.scale(m);
        grav.div(dens);

        acc.add(grav);

        vel.add(acc.scaleDup(dt * m).divDup(dens));
        pos.add(vel.scaleDup(dt));
        // pos.clampVec(bounds.min, bounds.max);

        inline for (0..3) |a| {
            const bmin = bounds.min.v[a];
            const bmax = bounds.max.v[a];

            if (pos.v[a] < bmin) {
                pos.v[a] = bmin;
                vel.v[a] *= boundDamping;
            }

            if (pos.v[a] > bmax) {
                pos.v[a] = bmax;
                vel.v[a] *= boundDamping;
            }
        }
        //pred.* = pos.*;
    }
}

fn predict(particles: ParticleList, dt: f32) void {
    var position = particles.items(.position);
    var predicted = particles.items(.predictedPos);
    var velocity = particles.items(.velocity);
    var acceleration = particles.items(.acceleration);
    var mass = particles.items(.mass);
    var density = particles.items(.density);

    for (
        0..particles.len,
        acceleration,
        velocity,
        position,
        predicted,
        mass,
        density,
    ) |p, *acc, *vel, pos, *pred, m, dens| {
        _ = acc;
        _ = p;

        var grav = particleWorld.gravity;
        grav.scale(m * dens);
        vel.add(grav.scaleDup(dt * m).divDup(dens));
        pred.* = pos.addDup(vel.scaleDup(dt));
    }
}

fn collideBounds(pos: *Vec4f, vel: *Vec4f) void {
    _ = vel;
    _ = pos;
    const bounds = particleWorld.bounds;
    _ = bounds;
}

fn updateDensity(particles: ParticleList, pw: ParticleWorld) void {
    var position = particles.items(.predictedPos);
    var density = particles.items(.density);
    var pressure = particles.items(.pressure);
    var mass = particles.items(.mass);
    var color = particles.items(.color);

    for (0..particles.len, density, position, color, pressure) |a, *adensity, apos, *acolor, *apressure| {
        adensity.* = 0;
        for (0..particles.len, position, mass) |b, bpos, bmass| {
            if (a == b)
                continue;

            const dir = apos.subDup(bpos);
            const distsq = dir.length3Sqr();

            if (distsq < smoothingRadiusSqr) {
                adensity.* += bmass * pw.POLY6 * std.math.pow(f32, smoothingRadiusSqr - distsq, 3.0);
            }
        }

        apressure.* = pw.gas * (adensity.* - pw.restingDensity);

        acolor.* = Vec4f.init(
            @mod(adensity.*, 0.75),
            @mod(adensity.*, 0.5),
            0.5 + @mod(adensity.*, 0.25),
            1.0,
        );
    }
}

fn updatePressureForces(particles: ParticleList, pw: ParticleWorld) void {
    var position = particles.items(.predictedPos);
    var density = particles.items(.density);
    var pressure = particles.items(.pressure);
    var velocity = particles.items(.velocity);
    var acceleration = particles.items(.acceleration);
    const bounds = pw.bounds;

    for (0..particles.len, position, pressure, velocity, acceleration, density) |a, apos, apressure, avelocity, *aaccel, adensity| {
        var pressAcc = Vec4f.zero();
        var viscAcc = Vec4f.zero();

        for (0..particles.len, position, pressure, density, velocity) |b, bpos, bpressure, bdensity, bvelocity| {
            if (a == b)
                continue;

            const dir = apos.subDup(bpos);
            const normDir = dir.normalized3();
            const dist = dir.length3();

            if (dist < smoothingRadius) {
                const radiusDiff = smoothingRadius - dist;
                const dv = bvelocity.subDup(avelocity);
                const pressScalar = (apressure + bpressure) / (2.0 * bdensity) * pw.SPIKY_GRAD * std.math.pow(f32, radiusDiff, 3.0);

                pressAcc.add(normDir.scaleDup(-pressScalar));

                viscAcc.add(dv.scaleDup(pw.VISC).divDup(bdensity * pw.VISC_LAP * radiusDiff));
            }
        }

        // calculate boundary forces
        inline for (0..3) |axis| {
            const bmin = bounds.min.v[axis];
            const bmax = bounds.max.v[axis];

            const dist = [_]f32{ apos.v[axis] - bmin, bmax - apos.v[axis] };

            inline for (0..1) |b| {
                if (dist[b] < smoothingRadius) {
                    const radiusDiff = smoothingRadius - dist[b];
                    const dv = avelocity;
                    const pressScalar = (apressure) / (2.0 * adensity) * pw.SPIKY_GRAD * std.math.pow(f32, radiusDiff, 3.0);

                    pressAcc.v[axis] += (@as(f32, @floatFromInt(b)) * 1.0 - 1.0) * -pressScalar;

                    viscAcc.v[axis] += (dv.v[axis] * pw.VISC) / (adensity * pw.VISC_LAP * radiusDiff);
                }
            }
        }

        aaccel.* = pressAcc;
        aaccel.add(viscAcc);
    }
}

var lightPos = Vec4f.init(10, -1, 4, 0);
const MaxSteps: usize = 50;
const maxDist: f32 = 1000.0;
const minSurfaceDist = 0.001;

pub fn renderPixel(view: *const Mat44f, fragCoord: Vec4f, viewport: Vec4f) Color {
    //_=view;
    //var rayOrigin = Vec4f.init(0,1,0,0);
    var rayOrigin = view.position();
    _ = rayOrigin;

    // // center in viewport
    const ray = Vec4f.init(
        (fragCoord.x() - viewport.x() * 0.5) / viewport.y(),
        (viewport.y() * 0.5 - fragCoord.y()) / viewport.y(),
        1,
        0,
    ).normalized();
    _ = ray;

    // const dist = march(rayOrigin, ray, MaxSteps);

    // if (dist >= maxDist)
    return Color.fromNormalVec4f(Vec4f.one().scaleDup(0.25));
    // const p = rayOrigin.addDup(ray.scaleDup(dist));
    // const n = sdf.normal2(p, dist, map);

    // const c = lightAt(p, n, lightPos);

    // const vc = Vec4f.splat(c);

    // return Color.fromNormalVec4f(uncharted2_filmic(vc));
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
    _ = m;
    _ = mv;
    _ = pixelLocal;
    _ = material;
    _ = uv;
    _ = normal;
    // var c = color.addDup(
    //   engine.Vec4f.init(
    //     (std.math.sin(uv.x*uv.y*1000)+1/2),
    //     (std.math.cos(uv.y*1000)+1/2),
    //     0,1)
    //   );

    var c = color;
    _ = pixel;
    _ = mvp;
    //var c = material.texture.sample(uv.x(), uv.y());
    //var c = material.texture.sampleBilinear(uv.x(), uv.y());

    //const l = @max(normal.dot3(material.lightDirection) * material.lightIntensity, 0.4);

    //c.scale(l);

    return c;
    //return uncharted2_filmic(c);
    //kreturn reinhard(c);
    //return c.scaleDup(l);
}
