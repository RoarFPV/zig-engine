const std = @import("std");
const math = std.math;
const Vec4f = @import("vector.zig").Vec4f;
const lerp = @import("interp.zig").lerp;

// https://www.shadertoy.com/view/Xds3zN

pub fn Sphere(point: Vec4f, radius: f32) f32 {
    return point.length() - radius;
}

pub fn Box(point: Vec4f, box: Vec4f) f32 {
    const q = point.abs().subDup(box);
    return q.maxScalar(0.0).length() + @min(q.maxElement(), 0.0);
}

pub fn Cube(point: Vec4f, size: f32) f32 {
    return Box(point, Vec4f.splat(size));
}

pub fn Plane(point: Vec4f, n: Vec4f, h: f32) f32 {
    return point.dot(n) + h;
}
const WorldDistanceFunc = *const fn (p: Vec4f) f32;

pub fn normal(p: Vec4f, d: f32, map: WorldDistanceFunc) Vec4f {
    const n = Vec4f.init(d - map(p.subDup(Vec4f.init(0.01, 0, 0, 0))), d - map(p.subDup(Vec4f.init(0, 0.01, 0, 0))), d - map(p.subDup(Vec4f.init(0, 0, 0.01, 0))), 0.0);

    return n.normalized();
}

pub fn normal2(p: Vec4f, d: f32, map: WorldDistanceFunc) Vec4f {
    _ = d;
    // inspired by tdhooper and klems - a way to prevent the compiler from inlining map() 4 times
    var n = Vec4f.splat(0.0);
    for (0..4) |i| {
        var e = Vec4f.init(@as(f32, @floatFromInt(((i + 3) >> 1) & 1)), @as(f32, @floatFromInt((i >> 1) & 1)), @as(f32, @floatFromInt(i & 1)), 0);
        e.scale(2);
        e.subScalar(1.0);
        e.scale(0.5773);

        e.scale(map(p.addDup(e.scaleDup(0.0005))));

        n.add(e);
        //if( n.x+n.y+n.z>100.0 ) break;
    }
    n.normalize();
    return n;
}

pub const Ops = struct {
    pub inline fn min(a: f32, b: f32) f32 {
        return @min(a, b);
    }

    pub inline fn sub(a: f32, b: f32) f32 {
        return @max(-a, b);
    }

    pub inline fn intersect(a: f32, b: f32) f32 {
        return @max(a, b);
    }

    pub const Smooth = struct {
        pub inline fn min(a: f32, b: f32, k: f32) f32 {
            const h = math.clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
            return lerp(b, a, h) - k * h * (1.0 - h);
        }

        pub inline fn sub(a: f32, b: f32, k: f32) f32 {
            const h = math.clamp(0.5 - 0.5 * (b + a) / k, 0.0, 1.0);
            return lerp(b, -a, h) + k * h * (1.0 - h);
        }

        pub inline fn intersect(a: f32, b: f32, k: f32) f32 {
            const h = math.clamp(0.5 - 0.5 * (b - a) / k, 0.0, 1.0);
            return lerp(b, a, h) + k * h * (1.0 - h);
        }
    };
};
// float sdBox( vec3 p, vec3 b )
// {
//   vec3 q = abs(p) - b;
//   return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
// }
