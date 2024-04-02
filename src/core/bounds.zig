// const std = @import("std");
const Vec4f = @import("vector.zig").Vec4f;

pub const Bounds = struct {
    min: Vec4f,
    max: Vec4f,

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
        self.min.sub(Vec4f.one().divDup(2.0));
        self.min.ceil();

        self.max.sub(Vec4f.one().divDup(2.0));
        self.max.ceil();
    }

    pub fn size(self: Bounds) Vec4f {
        return self.max.subDup(self.min);
    }
};
