// const std = @import("std");
const Vec4f = @import("vector.zig").Vec4f;
const std = @import("std");

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
        const color = Color{ .color = [4]u8{ 255, 255, 255, 255 } };
        return color;
    }

    pub fn black() Color {
        const color = Color{ .color = [4]u8{ 0, 0, 0, 255 } };
        return color;
    }

    pub fn red() Color {
        return .{ .color = [4]u8{ 255, 0, 0, 255 } };
    }

    pub fn green() Color {
        return .{ .color = [4]u8{ 0, 255, 0, 255 } };
    }

    pub fn blue() Color {
        return .{ .color = [4]u8{ 0, 0, 255, 255 } };
    }

    pub fn init(cr: u8, cg: u8, cb: u8, ca: u8) Color {
        return Color{ .color = [4]u8{ cr, cg, cb, ca } };
    }

    pub fn fromNormal(cr: f32, cg: f32, cb: f32, ca: f32) Color {
        return Color.init(
            @as(u8, @intFromFloat(std.math.clamp(cr, 0.0, 1.0) * 255 + 0.5)),
            @as(u8, @intFromFloat(std.math.clamp(cg, 0.0, 1.0) * 255 + 0.5)),
            @as(u8, @intFromFloat(std.math.clamp(cb, 0.0, 1.0) * 255 + 0.5)),
            @as(u8, @intFromFloat(std.math.clamp(ca, 0.0, 1.0) * 255 + 0.5)),
        );
    }

    pub fn fromNormalVec4f(vec: Vec4f) Color {
        const vc = vec.clamped01().scaleDup(255);
        //const xc: @Vector(4, i32) = @intFromFloat(vc.v);
        return .{ .color = [4]u8{
            @intFromFloat(vc.v[0]),
            @intFromFloat(vc.v[1]),
            @intFromFloat(vc.v[2]),
            @intFromFloat(vc.v[3]),
        } };
    }

    pub fn toNormalVec4f(c: Color) Vec4f {
        return Vec4f.init(
            @as(f32, @floatFromInt(c.r())) / 255.0,
            @as(f32, @floatFromInt(c.g())) / 255.0,
            @as(f32, @floatFromInt(c.b())) / 255.0,
            @as(f32, @floatFromInt(c.a())) / 255.0,
        );
    }
};
