// compile in ubuntu:
// $ zig build-exe paint.zig --library SDL2 --library SDL2main --library c -isystem "/usr/include" --library-path "/usr/lib/x86_64-linux-gnu"

const std = @import("std");
const warn = std.debug.print;
const fmt = std.fmt;
const assert = @import("std").debug.assert;
const math = std.math;

const Vec4f = @import("../core/vector.zig").Vec4f;
const Mat44f = @import("../core/matrix.zig").Mat44f;
const Color = @import("../core/color.zig").Color;

const Profile = @import("../core/profiler.zig").Profile;
pub const Font = @import("font.zig").Font;

const Mesh = @import("mesh.zig").Mesh;
const pixelbuffer = @import("pixel_buffer.zig");
const PixelBuffer = pixelbuffer.PixelBuffer;
const PixelRenderer = pixelbuffer.PixelRenderer;
const sdf = @import("../core/sdf.zig");

pub const material = @import("material.zig");
pub const Material = material.Material;

const tracy = @import("../tracy.zig");
const trace = tracy.trace;

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

var stats = Stats.init();
var viewport = Vec4f.zero();
var renderBounds: Bounds = undefined;

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

    colorRenderer = PixelRenderer(Color).init(&colorBuffer[currentBuffer]);

    renderBounds = Bounds.init(Vec4f.zero(), viewport);
    viewport = Vec4f.init(@as(f32, @floatFromInt(colorBuffer[0].w)), @as(f32, @floatFromInt(colorBuffer[0].h)), 0, 0);
}

pub fn shutdown() void {
    colorBuffer[0].deinit();
    colorBuffer[1].deinit();
    depthBuffer[0].deinit();
    depthBuffer[1].deinit();
}

pub fn beginFrame(profiler: ?*Profile) *u8 {
    const zone = trace(@src());
    defer zone.end();

    profile = profiler;
    var pprof = profile.?.beginSample("render.beginFrame", 0);
    defer profile.?.endSample(pprof);

    stats.reset();

    currentBuffer = ~currentBuffer;
    colorRenderer.setBuffer(&colorBuffer[currentBuffer]);

    return &colorBuffer[currentBuffer].bufferStart().color[0];
}

pub fn endFrame() void {}

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
    const px = shader.projectionShader(
        proj,
        shader.vertexShader(mv, 0, point, shader),
        viewport,
        shader,
    );

    //   shader.vertexShader(mvp, 0, point, shader);

    if (px.x() >= 0 and px.x() <= 1000 and px.y() >= 0 and px.y() <= 1000) {
        const pc = color;

        writePixel(
            @as(i32, @intFromFloat(px.x())),
            @as(i32, @intFromFloat(px.y())),
            px.z(),
            Color.fromNormalVec4f(pc.clamped01()),
        );
    }
}

///
pub fn drawWorldLine(proj: *const Mat44f, mv: *const Mat44f, start: Vec4f, end: Vec4f, color: Vec4f, shader: *Material) void {
    const spx = shader.projectionShader(proj, shader.vertexShader(mv, 0, start, shader), viewport, shader);
    const epx = shader.projectionShader(proj, shader.vertexShader(mv, 0, end, shader), viewport, shader);
    const pc = color;

    if (spx.x() >= 0 and spx.x() <= 1000 and spx.y() >= 0 and spx.y() <= 1000 and
        epx.x() >= 0 and epx.x() <= 1000 and epx.y() >= 0 and epx.y() <= 1000)
    {
        drawLine(
            @as(i32, @intFromFloat(spx.x())),
            @as(i32, @intFromFloat(spx.y())),
            @as(i32, @intFromFloat(epx.x())),
            @as(i32, @intFromFloat(epx.y())),
            Color.fromNormalVec4f(pc),
        );
    }
}

pub fn drawProgress(x: i16, y: i16, max_width: f32, value: f32, max_value: f32) void {
    const cs = std.math.clamp(value, 0.0, max_value) / max_value;
    // const cs2 = cs*cs;
    drawLine(
        x,
        y,
        @as(c_int, @intFromFloat(cs * max_width)),
        y,
        Color.fromNormal(cs, (1 - cs), 0.2, 1),
    );
}

pub fn writePixelNormal(x: i32, y: i32, z: f32, c: Vec4f) void {
    colorBuffer[currentBuffer].write(x, y, Color.fromNormal(c.x(), c.y(), c.z(), c.w()));
    depthBuffer[currentBuffer].write(x, y, z);
}

pub fn writePixel(x: i32, y: i32, z: f32, c: Color) void {
    colorBuffer[currentBuffer].write(x, y, c);
    depthBuffer[currentBuffer].write(x, y, z);
}

const RenderPixelFunc = *const fn (view: *const Mat44f, fragCoord: Vec4f, viewport: Vec4f) Color;

pub fn renderFrame(view: *const Mat44f, renderPixel: RenderPixelFunc) void {
    const w = @as(usize, @intCast(colorBuffer[currentBuffer].w));
    const h = @as(usize, @intCast(colorBuffer[currentBuffer].h));

    for (0..h) |uy| {
        for (0..w) |ux| {
            const x = @as(i32, @intCast(ux));
            const y = @as(i32, @intCast(uy));

            //const index = colorBuffer.pxIndex(x,y);
            const fragCoord = Vec4f.init(@as(f32, @floatFromInt(x)), @as(f32, @floatFromInt(y)), 0.0, 0.0);

            const color = renderPixel(view, fragCoord, viewport);

            writePixel(x, y, 1, color);
        }
    }
}
