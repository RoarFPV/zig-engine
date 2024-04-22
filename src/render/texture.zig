const std = @import("std");
const warn = std.debug.print;
const fmt = std.fmt;
const assert = @import("std").debug.assert;
const math = std.math;
const io = std.io;

const Vec4f = @import("../core/vector.zig").Vec4f;
const Mat44f = @import("../core/matrix.zig").Mat44f;
const Profile = @import("../core/profiler.zig").Profile;

pub const trace = @import("../tracy.zig").trace;

pub const Format = enum {
    GRAY8,
    RGB8,
    RGBA8,
    F32,
};

pub const Texture = struct {
    width: f32,
    height: f32,
    iwidth: i32,
    iheight: i32,
    colors: []Vec4f,
    format: Format,
    widthBytes: u32,
    heightBytes: u32,
    pixelWidth: u32,

    pub const Error = error{InvalidImageFormat};

    inline fn normalizeColor(color: u8) f32 {
        return @as(f32, @floatFromInt(color)) / 255;
    }

    pub fn init(allocator: *std.mem.Allocator, f: Format, w: u32, h: u32, pixelBytes: u8, data: []u8) !Texture {
        var colors = try allocator.alloc(Vec4f, w * h);

        for (0..(w * h)) |s| {
            const ps = (s * pixelBytes);
            const source = data[ps..(ps + pixelBytes)];

            colors[s] = switch (f) {
                .GRAY8 => Vec4f.init(normalizeColor(source[0]), normalizeColor(source[0]), normalizeColor(source[0]), 1),
                .RGB8 => Vec4f.init(normalizeColor(source[0]), normalizeColor(source[1]), normalizeColor(source[2]), 1),
                .RGBA8 => Vec4f.init(
                    normalizeColor(source[0]),
                    normalizeColor(source[1]),
                    normalizeColor(source[2]),
                    normalizeColor(source[3]),
                ),
                .F32 => return error.InvalidImageFormat,
            };
        }

        return Texture{
            .format = f,
            .width = @as(f32, @floatFromInt(w)),
            .height = @as(f32, @floatFromInt(h)),
            .iwidth = @as(i32, @intCast(w)),
            .iheight = @as(i32, @intCast(h)),
            .colors = colors,
            .widthBytes = w,
            .heightBytes = h,
            .pixelWidth = pixelBytes,
        };
    }

    pub fn initVec(allocator: *std.mem.Allocator, f: Format, w: u32, h: u32, data: []Vec4f) !Texture {
        _ = allocator;
        return Texture{
            .format = f,
            .width = @as(f32, @floatFromInt(w)),
            .height = @as(f32, @floatFromInt(h)),
            .iwidth = @as(i32, @intCast(w)),
            .iheight = @as(i32, @intCast(h)),
            .colors = data,
            .widthBytes = w,
            .heightBytes = h,
            .pixelWidth = @sizeOf(Vec4f),
        };
    }

    pub fn sample(self: Texture, x: f32, y: f32) Vec4f {
        const tx = @as(usize, @intFromFloat(std.math.clamp(x, 0, 1) * (self.width - 1)));
        const ty = @as(usize, @intFromFloat(std.math.clamp(y, 0, 1) * (self.height - 1)));
        const index = ty * self.widthBytes + tx;

        return self.colors[index];
    }

    pub fn sampleIndex(self: Texture, index: usize) Vec4f {
        return self.colors[index];
    }

    pub fn samplePixel(self: Texture, x: i32, y: i32) Vec4f {
        const tx = @as(u32, @intCast(std.math.clamp(x, 0, self.iwidth - 1)));
        const ty = @as(u32, @intCast(std.math.clamp(y, 0, self.iheight - 1)));
        const index = ty * self.widthBytes + tx;

        return self.colors[index];
    }

    pub fn sampleBilinear(self: Texture, x: f32, y: f32) Vec4f {
        const zone = trace(@src());
        defer zone.end();

        const tx = @as(i32, @intFromFloat(std.math.clamp(x, 0, 1) * (self.width - 1)));
        const ty = @as(i32, @intFromFloat(std.math.clamp(y, 0, 1) * (self.height - 1)));

        const dx = x * self.width - @as(f32, @floatFromInt(tx));
        const dy = y * self.height - @as(f32, @floatFromInt(ty));

        const s0 = self.samplePixel(tx, ty);
        const s1 = self.samplePixel(tx + 1, ty);
        const s2 = self.samplePixel(tx, ty + 1);
        const s3 = self.samplePixel(tx + 1, ty + 1);

        return Vec4f.blerpDup(s0, s1, s2, s3, dx, dy);
    }

    pub fn sampleBiCubic(self: Texture, x: f32, y: f32) Vec4f {
        return self.sample(x, y);
    }
};
