// platform imports
const std = @import("std");
const fmt = std.fmt;
const warn = std.debug.print;
const assert = std.debug.assert;
const Timer = std.time.Timer;

// engine imports
pub const sys = @import("system/sys_sdl.zig");
//pub const render = @import("render/raster_cpu_jobs.zig");
// pub const render = @import("render/raytracer_cpu.zig");
pub const render = @import("render/raster_cpu_spans.zig");
pub const input = @import("system/sys_input.zig");
pub const game = @import("game/game_cubes.zig");
// pub const game = @import("game/game_voxel.zig");
// pub const game = @import("game/game_particle.zig");
//pub const game = @import("game/game_oblivious.zig");

pub const Mat44f = @import("core/matrix.zig").Mat44f;
pub const Vec4f = @import("core/vector.zig").Vec4f;
pub const Profile = @import("core/profiler.zig").Profile;
pub const Sampler = @import("core/profiler.zig").Sampler;
pub const Color = @import("core/color.zig").Color;
pub const blend = @import("core/interp.zig");

pub const Mesh = @import("render/mesh.zig").Mesh;
pub const Terrain = @import("render/terrain.zig").Terrain;

pub const tracy = @import("tracy.zig").tracy;
pub const trace = tracy.trace;

const tools = @import("tools.zig");

pub var stdout = std.io.getStdOut();
pub var deltaTimeNs: u64 = 0;
pub var deltaTimeSec: f32 = 0;

var bufferAllocator = std.heap.page_allocator;
var textureAllocator = std.heap.page_allocator;
pub var font: render.Font = undefined;

var render3d: bool = true;
var renderSingleFrame: bool = false;

// pub const systemConfig = sys.Config{
//     .windowWidth = 1024,
//     .windowHeight = 768,
//     .renderWidth = 1024,
//     .renderHeight = 768,
//     .maxFps = 60,
//     .fullscreen = false,
// };

// pub const systemConfig = sys.Config{
//     .windowWidth = 1920,
//     .windowHeight = 1080,
//     .renderWidth = 1920,
//     .renderHeight = 1080,
//     .maxFps = 60,
//     .fullscreen = false,
// };

pub const systemConfig = sys.Config{
    //.windowWidth = 1024,
    .windowWidth = 1260,
    .windowHeight = 768,
    // .renderWidth = 328,
    // .renderHeight = 200,
    // .renderWidth = 640,
    .renderWidth = 787,
    .renderHeight = 480,
    .maxFps = 60,
    .fullscreen = false,
};

// pub const systemConfig = sys.Config{
//     .windowWidth = 1260,
//     .windowHeight = 768,
//     .renderWidth = 1260,
//     .renderHeight = 768,
//     .maxFps = 60,
//     .fullscreen = false,
// };

var profileId: u1 = 0;

var profiles: [2]Profile = undefined;

fn createAsciiTable(comptime s: usize, comptime n: usize) [n]u8 {
    assert(s < 255);
    assert(n < 255);

    var result: [n]u8 = undefined;
    for (s..n) |i| {
        const c = @as(u8, i);
        result[i - s] = c;
        //@compileLog("ascii[{d}]={d}", i-s, c);
    }
    return result;
}

const asciiTable = createAsciiTable(16 * 2, 32);

pub fn swapProfile() *Profile {
    profileId = ~profileId;
    return currentProfile();
}

pub fn currentProfile() *Profile {
    return &profiles[profileId];
}

pub fn nextProfile() *Profile {
    return &profiles[~profileId];
}

pub fn main() !void {
    profiles = [2]Profile{ try Profile.init(&bufferAllocator), try Profile.init(&bufferAllocator) };

    var profiler = currentProfile();
    profiler.nextFrame();

    var lastProfile = nextProfile();
    lastProfile.nextFrame();

    try sys.init(systemConfig);
    defer sys.shutdown();

    try render.init(systemConfig.renderWidth, systemConfig.renderHeight, &bufferAllocator, profiler);
    defer render.shutdown();

    const fontTex = try tools.TgaTexLoader.importTGAFile(&textureAllocator, "../../assets/mbf_small_7x7.tga");

    font = render.Font{
        .glyphWidth = 7,
        .glyphHeight = 7,
        .texture = fontTex,
    };

    //defer textureAllocator.reset();

    const bufferLineSize = render.bufferLineSize();

    var quit = false;
    const targetFrameTimeNs = @as(f32, @floatFromInt(sys.targetFrameTimeMs() * 1_000_000));

    _ = try game.init();

    // var mainSampler:Sampler = undefined;
    // var systemSampler:Sampler = undefined;
    // var gameSampler:Sampler = undefined;
    // var renderSampler:Sampler = undefined;

    //std.debug.print("{s}\n", .{&asciiTable});

    while (!quit) {
        const mainTrace = trace(@src());
        defer mainTrace.end();

        const startTime = try std.time.Instant.now();

        {
            tracy.frameMark();
            var el = Sampler.initAndBegin(profiler, "engine.main", 0);
            defer el.end();

            {
                var supdate = Sampler.initAndBegin(profiler, "system.update", 1);
                defer supdate.end();

                quit = !sys.beginUpdate();
            }

            const b = render.beginFrame(profiler);
            {
                var c = Sampler.initAndBegin(profiler, "game.update", 2);
                defer c.end();
                if (!game.update())
                    break;
            }

            {
                var ef = Sampler.initAndBegin(profiler, "render.endframe", 3);
                defer ef.end();
                render.endFrame();
            }

            {
                var srt = Sampler.initAndBegin(profiler, "engine.render.profile", 0);
                defer srt.end();

                if (lastProfile.hasSamples()) {
                    // render.drawString(&font, "0123456789\n", 10, 50, Vec4f.one());
                    displayProfileUi(
                        lastProfile,
                        2,
                        2,
                        5,
                        systemConfig.renderWidth - 5,
                        8,
                        targetFrameTimeNs,
                    );
                    //render.drawProgress(2, 2, 100, @intToFloat(f32, lastProfile.sampleTime(1)), targetFrameTimeNs );
                    //render.drawProgress(2, 5, 100, @intToFloat(f32, lastProfile.sampleTime(2)), targetFrameTimeNs );
                }
            }

            {
                var srt = Sampler.initAndBegin(profiler, "system.render.present", 2);
                defer srt.end();

                sys.updateRenderTexture(b, bufferLineSize);
                sys.renderPresent();
            }

            {
                var sup = Sampler.initAndBegin(profiler, "system.wait", 1);
                defer sup.end();

                _ = sys.endUpdate();
            }
        }

        lastProfile = profiler;
        profiler = swapProfile();
        profiler.nextFrame();

        const endTime = try std.time.Instant.now();

        deltaTimeNs = endTime.since(startTime);

        deltaTimeSec = @floatFromInt(deltaTimeNs);
        deltaTimeSec /= (1000.0 * 1000.0 * 1000.0);
    }

    // try currentProfile().jsonFileWrite(bufferAllocator, "prof.json");

    // bufferAllocator.reset();
}

const colors = [_]Color{
    Color.white(),
    Color.red(),
    Color.green(),
    Color.blue(),
    Color.fromNormal(0.5, 0.2, 1, 1),
};

pub fn displayProfileUi(self: *Profile, x: i32, y: i32, lineSize: i32, maxWidth: f32, maxDepth: u8, targetNs: f32) void {
    const mainSample = self.samples.items[1];
    // const totalBegin = @intToFloat(f32, mainSample.begin - self.frameStartTime);
    const totalEnd: f32 = @floatFromInt(mainSample.end.since(self.frameStartTime));
    const totalTime = targetNs * 2; //@intToFloat(f32, mainSample.end-mainSample.begin);

    const ny = y + font.glyphHeight + lineSize;
    // y += 10;

    const targetx = x + @as(i32, @intFromFloat((targetNs / totalTime) * maxWidth));
    render.drawLine(
        targetx,
        ny - 2,
        targetx,
        ny + lineSize * @as(i32, @intCast(maxDepth)),
        Color.fromNormal(0.5, 0.0, 0.2, 0.7),
    );

    const targetmidx = x + @as(i32, @intFromFloat(((targetNs / 2) / totalTime) * maxWidth));
    render.drawLine(
        targetmidx,
        ny - 2,
        targetmidx,
        ny + lineSize * @as(i32, @intCast(maxDepth)),
        Color.fromNormal(0.3, 0.3, 0.3, 0.7),
    );

    for (self.samples.items, 0..) |sample, i| {
        if (i > self.nextSample)
            break;

        if (sample.depth >= maxDepth)
            continue;

        if (i == 0 or sample.begin.since(self.frameStartTime) <= 0 or sample.begin.order(self.frameStartTime) == std.math.Order.lt)
            continue;

        const begin: f32 = @floatFromInt(sample.begin.since(self.frameStartTime));
        const end: f32 = @floatFromInt(sample.end.since(self.frameStartTime));
        // const duration = end-begin;

        // const cs = std.math.clamp(blend.invLerp(f32, 0, targetNs, totalEnd), 0.0, 1.0);

        const startx = x + @as(i32, @intFromFloat(std.math.clamp(begin / totalTime, 0.0, 4.0) * maxWidth));
        const finishx = x + @as(i32, @intFromFloat(std.math.clamp(end / totalTime, 0.0, 4.0) * maxWidth));
        const ystart = ny + lineSize * @as(i32, @intCast(sample.depth));

        // std.debug.warn("{} x: {}, fx:{}, y:{}, b:{}, e:{}, d:{}, t:{}\n", .{sample.depth, startx, finishx, ystart, begin, end, duration, sample.tag});

        render.drawLine(startx, ystart, finishx, ystart, colors[sample.color]);
    }

    var buf: [64]u8 = undefined;
    const b = std.fmt.bufPrint(&buf, "{d:.3}", .{totalEnd / 1_000_000.0}) catch return;
    render.drawString(&font, b, x, y, Color.white().toNormalVec4f());
}

// 0 x: 2, fx:202, y:2, b:2.61e+02, e:7.813682e+06, d:7.813421e+06, t:engine.main
//   1 x: 2, fx:2, y:5, b:3.85e+02, e:6.299e+03, d:5.914e+03, t:system.update
//   1 x: 2, fx:66, y:5, b:6.457e+03, e:2.516793e+06, d:2.510336e+06, t:render.beginFrame
//   1 x: 66, fx:181, y:5, b:2.517428e+06, e:7.022163e+06, d:4.504735e+06, t:game.update
//     2 x: 66, fx:180, y:8, b:2.526334e+06, e:6.978469e+06, d:4.452135e+06, t:render.mesh.draw
//   1 x: 181, fx:190, y:5, b:7.022379e+06, e:7.379825e+06, d:3.57446e+05, t:engine.render.profile
//   1 x: 190, fx:201, y:5, b:7.380564e+06, e:7.810696e+06, d:4.30132e+05, t:system.render.present
//   1 x: 201, fx:201, y:5, b:7.811119e+06, e:7.81359e+06, d:2.471e+03, t:system.wait
