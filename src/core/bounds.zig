// const std = @import("std");
const Vec4f = @import("vector.zig").Vec4f;
const std = @import("std");

pub const Bounds = struct {
    min: Vec4f,
    max: Vec4f,

    pub fn initInfinity() Bounds {
        return Bounds{
            .min = Vec4f.splat(std.math.inf(f32)),
            .max = Vec4f.splat(-std.math.inf(f32)),
        };
    }

    pub fn init(min: Vec4f, max: Vec4f) Bounds {
        return Bounds{
            .min = min,
            .max = max,
        };
    }

    pub fn add(self: *Bounds, point: Vec4f) void {
        self.min = self.min.min(point);
        self.max = self.max.max(point);
    }

    pub fn limit(self: *Bounds, l: Bounds) void {
        self.min = Vec4f.max(l.min, self.min);
        self.max = Vec4f.min(l.max, self.max);
    }

    pub fn topLeftHandLimit(self: *Bounds) void {
        const half = Vec4f.splat(0.5);
        self.min.sub(half);
        self.min.ceil();

        self.max.sub(half);
        self.max.ceil();
    }

    pub fn size(self: Bounds) Vec4f {
        return self.max.subDup(self.min);
    }

    pub fn halfSize(self: Bounds) Vec4f {
        return self.size().divDup(2);
    }

    pub fn center(self: Bounds) Vec4f {
        return self.max
            .subDup(self.min)
            .divDup(2)
            .addDup(self.min);
    }
};
