const std = @import("std");
const math = std.math;
const interp = @import("interp.zig");

pub const Vector3f = @Vector(3, f32);
const Vector2f = @Vector(2, f32);

pub fn Vec4(comptime eT: type) type {
    return struct {
        const Self = @This();
        pub const Vector4 = @Vector(4, eT);
        pub const ElementType = eT;

        v: Vector4,

        pub inline fn init(_x: eT, _y: eT, _z: eT, _w: eT) Self {
            return .{ .v = .{ _x, _y, _z, _w } };
        }

        pub inline fn initXY(sx: eT, sy: eT) Self {
            return .{ .v = .{ sx, sy, 0, 1 } };
        }

        const Axis = enum(u8) {
            X = 0,
            Y = 1,
            Z = 2,
            W = 3,
        };

        pub inline fn initVector(v: Vector4) Self {
            return .{ .v = v };
        }

        pub inline fn splat(s: eT) Self {
            return .{ .v = .{ s, s, s, s } };
        }

        pub inline fn splat3(s: eT, _w: eT) Self {
            return .{ .v = .{ s, s, s, _w } };
        }

        pub fn zero() Self {
            return Self.init(0, 0, 0, 0);
        }

        pub fn one() Self {
            return Self.init(1, 1, 1, 1);
        }

        pub fn forward() Self {
            return Self.init(0, 0, 1, 1);
        }

        pub fn up() Self {
            return Self.init(0, 1, 0, 1);
        }

        pub fn right() Self {
            return Self.init(1, 0, 0, 1);
        }

        // TODO: generate swizzle combinations?
        pub inline fn x(self: Self) eT {
            return self.v[@intFromEnum(Axis.X)];
        }

        pub inline fn y(self: Self) eT {
            return self.v[@intFromEnum(Axis.Y)];
        }

        pub inline fn z(self: Self) eT {
            return self.v[@intFromEnum(Axis.Z)];
        }

        pub inline fn w(self: Self) eT {
            return self.v[@intFromEnum(Axis.W)];
        }

        pub inline fn setValue(self: *Self, element: Axis, value: eT) void {
            self.v[@intFromEnum(element)] = value;
        }

        pub inline fn setX(self: *Self, value: eT) void {
            self.setValue(Axis.X, value);
        }

        pub inline fn setY(self: *Self, value: eT) void {
            self.setValue(Axis.Y, value);
        }

        pub inline fn setZ(self: *Self, value: eT) void {
            self.setValue(Axis.Z, value);
        }

        pub inline fn setW(self: *Self, value: eT) void {
            self.setValue(Axis.W, value);
        }

        pub inline fn setXY(self: *Self, sx: eT, sy: eT) void {
            self.setValue(Axis.X, sx);
            self.setValue(Axis.Y, sy);
        }

        pub inline fn setXYZ(self: *Self, sx: eT, sy: eT, sz: eT) void {
            self.setValue(Axis.X, sx);
            self.setValue(Axis.Y, sy);
            self.setValue(Axis.Z, sz);
        }

        pub inline fn fromXY(self: *Self, sx: eT, sy: eT) void {
            self.v = .{ sx, sy, 0, 1 };
        }

        pub inline fn fromXYZ(self: *Self, sx: eT, sy: eT, sz: eT) void {
            self.v = .{ sx, sy, sz, 0 };
        }

        pub inline fn set(self: *Self, other: Self) void {
            self.v = other.v;
        }

        pub inline fn add(self: *Self, other: Self) void {
            self.v += other.v;
        }

        pub inline fn addScalar(self: *Self, scalar: eT) void {
            self.v += @as(Vector4, @splat(scalar));
        }

        pub inline fn addScalarDup(self: Self, scalar: eT) Self {
            return Self{ .v = self.v + @as(Vector4, @splat(scalar)) };
        }

        pub inline fn addDup(self: Self, other: Self) Self {
            return Self{ .v = self.v + other.v };
        }

        pub inline fn mul(self: *Self, other: Self) void {
            self.v *= other.v;
        }

        pub inline fn mulDup(self: Self, other: Self) Self {
            return Self{ .v = self.v * other.v };
        }

        pub inline fn scale(self: *Self, scalar: eT) void {
            self.v *= @as(Vector4, @splat(scalar));
        }

        pub fn eq(self: Self, other: Self, _epsilon: eT) bool {
            const vabs = @abs(self.v - other.v);
            return @reduce(.And, vabs < _epsilon);
        }

        pub fn equal(self: Self, other: Self) bool {
            return @reduce(.And, self.v == other.v);
        }

        pub inline fn min(lhs: Self, rhs: Self) Self {
            return .{ .v = @min(lhs.v, rhs.v) };
        }

        pub inline fn max(lhs: Self, rhs: Self) Self {
            return .{ .v = @max(lhs.v, rhs.v) };
        }

        pub inline fn maxElement(self: Self) eT {
            return @reduce(.Max, self.v);
        }

        pub inline fn minElement(self: Self) eT {
            return @reduce(.Min, self.v);
        }

        pub inline fn scaleDup(self: Self, scalar: eT) Self {
            return Self{ .v = self.v * @as(Vector4, @splat(scalar)) };
        }

        pub inline fn scale3Dup(self: Self, scalar: eT) Self {
            const vscalar = Vector4{ scalar, scalar, scalar, 1 };
            return .{ .v = self.v * vscalar };
        }

        pub inline fn div(self: *Self, scalar: eT) void {
            self.v /= @as(Vector4, @splat(scalar));
        }

        pub inline fn div3(self: *Self, scalar: eT) void {
            const vscalar = Vector4{ scalar, scalar, scalar, 1 };
            self.v /= vscalar;
        }

        pub inline fn divDup(self: Self, scalar: eT) Self {
            return Self{ .v = self.v / @as(Vector4, @splat(scalar)) };
        }

        pub inline fn divVec(self: *Self, other: Self) void {
            self.v /= other.v;
        }

        pub inline fn divVecDup(self: Self, other: Self) Self {
            return Self{ .v = self.v / other.v };
        }

        pub inline fn pow(self: *Self, other: Self) void {
            self.v = self.powDup(other).v;
        }

        pub inline fn powDup(self: *Self, other: Self) Self {
            return Self{ .v = .{
                std.math.pow(eT, self.v[0], other.v[0]),
                std.math.pow(eT, self.v[1], other.v[1]),
                std.math.pow(eT, self.v[2], other.v[2]),
                std.math.pow(eT, self.v[3], other.v[3]),
            } };
        }

        pub inline fn dot3(self: Self, other: Self) eT {
            const v: Vector4 = .{ 1, 1, 1, 0 };
            return @reduce(.Add, v * self.v * other.v);
            //return self.x() * other.x() + self.y() * other.y() + self.z() * other.z();
        }

        pub inline fn abs(self: Self) Self {
            return .{ .v = @abs(self.v) };
        }

        pub inline fn dot(self: Self, other: Self) eT {
            return @reduce(.Add, (self.v * other.v));
        }

        inline fn dotRaw(self: Self, vother: Vector4) eT {
            return @reduce(.Add, (self.v * vother));
        }

        pub inline fn cross3(self: Self, other: Self) Self {
            return Self.init(
                self.y() * other.z() - self.z() * other.y(),
                self.z() * other.x() - self.x() * other.z(),
                self.x() * other.y() - self.y() * other.x(),
                1,
            );
        }

        pub inline fn sub(self: *Self, other: Self) void {
            self.v -= other.v;
        }

        pub inline fn subDup(self: Self, other: Self) Self {
            return .{ .v = self.v - other.v };
        }

        pub inline fn subScalar(self: *Self, scalar: eT) void {
            self.v -= @as(Vector4, @splat(scalar));
        }

        pub inline fn neg(self: Self) Self {
            return .{ .v = -self.v };
        }
        /// Returns vector length
        pub inline fn length(self: Self) eT {
            return math.sqrt(self.dot(self));
        }

        /// return length of vector squared. Avoids `math.sqrt`
        pub inline fn lengthSqr(self: Self) eT {
            return self.dot(self);
        }

        /// make `length` of vector 1.0 while maintaining direction
        pub inline fn normalize(self: *Self) void {
            self.div(self.length());
        }

        /// constant version of normalize that returns a new `Self` with length of 1.0
        pub inline fn normalized(self: Self) Self {
            return self.divDup(self.length());
        }

        pub inline fn length3(self: Self) eT {
            return @sqrt(self.length3Sqr());
        }

        /// return length of vector squared. Avoids `math.sqrt`
        pub inline fn length3Sqr(self: Self) eT {
            return self.dot3(self);
        }

        /// make `length` of vector 1.0 while maintaining direction
        pub inline fn normalize3(self: *Self) void {
            self.div3(self.length3());
            // self.normaize();
        }

        /// constant version of normalize that returns a new `Self` with length of 1.0
        pub inline fn normalized3(self: Self) Self {
            const len = @as(Vector4, @splat(self.length3()));
            return .{ .v = self.v / len };
            // return Self.init(self.x() / len, self.y() / len, self.z() / len, self.w());
            //return self.normalized();
        }

        pub inline fn clamp01(self: *Self) void {
            self.clamp(0, 1);
        }

        pub inline fn clamped01(self: Self) Self {
            return .{ .v = self._clampedScalar(0, 1) };
        }

        inline fn _clampedScalar(self: Self, minVec4: eT, maxVec4: eT) Vector4 {
            const maxV = @as(Vector4, @splat(maxVec4));
            const minV = @as(Vector4, @splat(minVec4));

            return @min(maxV, @max(minV, self.v));
        }

        pub inline fn maxScalar(self: Self, scalar: eT) Self {
            const maxV = @as(Vector4, @as(Vector4, @splat(scalar)));
            return .{ .v = @max(self.v, maxV) };
        }

        pub inline fn clampScalar(self: *Self, minVec4: eT, maxVec4: eT) void {
            self.v = self._clampedScalar(minVec4, maxVec4);
        }

        inline fn _clampedVec(self: Self, minV: Self, maxV: Self) Vector4 {
            return @min(maxV.v, @max(minV.v, self.v));
        }

        pub inline fn clamped(self: Self, minVec4: Self, maxVec4: Self) Vec4f {
            return .{ .v = self._clampedVec(minVec4, maxVec4) };
        }

        pub inline fn clamp(self: *Self, minVec4: Self, maxVec4: Self) void {
            self.v = self._clampedVec(minVec4, maxVec4);
        }

        pub inline fn ceil(self: *Self) void {
            self.v = @ceil(self.v);
        }

        pub fn print(self: Self) void {
            std.debug.print(" [{} | {}]", .{ self.v, self.length() });
        }

        pub fn println(self: Self) void {
            self.print();
            std.debug.print("\n", .{});
        }

        pub fn triArea(a: Self, b: Self, c: Self) eT {
            return (c.x() - a.x()) * (b.y() - a.y()) -
                (c.y() - a.y()) * (b.x() - a.x());
        }

        pub fn triCoords(v0: Self, v1: Self, v2: Self, p: Self) Self {
            return Self.init(
                triArea(v1, v2, p),
                triArea(v2, v0, p),
                triArea(v0, v1, p),
                1,
            );
        }

        pub fn triBarycentericCoordsOld(v0: Self, v1: Self, v2: Self, p: Self) Self {
            return Self.init(
                triArea(v1, v2, p),
                triArea(v2, v0, p),
                triArea(v0, v1, p),
                1,
            );
        }

        pub fn triBarycentericCoords(a: Self, b: Self, c: Self, p: Self) Self {
            const v0 = b.subDup(a);
            const v1 = c.subDup(a);
            const v2 = p.subDup(a);
            const d00 = v0.dot(v0);
            const d01 = v0.dot(v1);
            const d11 = v1.dot(v1);
            const d20 = v2.dot(v0);
            const d21 = v2.dot(v1);
            const invDenom = 1 / (d00 * d11 - d01 * d01);

            const v = (d11 * d20 - d01 * d21) * invDenom;
            const _w = (d00 * d21 - d01 * d20) * invDenom;
            const u = 1 - v - _w;
            return Self.init(v, _w, u, 0);
        }

        pub fn triInterpStd(tri: Self, v0: Self, v1: Self, v2: Self, depth: eT, _w: eT) Self {
            return Self.init(
                (tri.x() * v0.x() + tri.y() * v1.x() + tri.z() * v2.x()) / depth,
                (tri.x() * v0.y() + tri.y() * v1.y() + tri.z() * v2.y()) / depth,
                (tri.x() * v0.z() + tri.y() * v1.z() + tri.z() * v2.z()) / depth,
                _w,
            );
        }

        pub fn triInterp(tri: Self, v0: Self, v1: Self, v2: Self, tw: f32) Self {
            const vx = Vector4{ v0.x(), v1.x(), v2.x(), tw };
            const vy = Vector4{ v0.y(), v1.y(), v2.y(), tw };
            const vz = Vector4{ v0.z(), v1.z(), v2.z(), tw };
            const vw = Vector4{ v0.w(), v1.w(), v2.w(), tw };

            return Self.init(
                tri.dotRaw(vx),
                tri.dotRaw(vy),
                tri.dotRaw(vz),
                tri.dotRaw(vw),
            );
        }

        pub fn triInterpArray(tri: Self, v: [3]Self, tw: f32) Self {
            return triInterp(
                tri,
                v[0],
                v[1],
                v[2],
                tw,
            );
        }

        pub fn triInterpArrayScale(tri: Self, v: [3]Self, tw: f32, s: f32) Self {
            return triInterp(
                tri,
                v[0],
                v[1],
                v[2],
                tw,
            ).scaleDup(s);
        }

        pub fn triInterpBaryArray(tri: Self, v: [4]Self) Self {
            return Self.init(
                tri.dotRaw(v[0].v),
                tri.dotRaw(v[1].v),
                tri.dotRaw(v[2].v),
                tri.dotRaw(v[3].v),
            );
        }

        pub fn lerpDup(from: Self, to: Self, d: eT) Self {
            return .{
                .v = interp.lerp(
                    Self.Vector4,
                    .{ 1, 1, 1, 1 },
                    from.v,
                    to.v,
                    Self.splat(d).v,
                ),
            };
        }

        pub fn blerpDup(fromx: Self, tox: Self, fromy: Self, toy: Self, dx: eT, dy: eT) Self {
            return lerpDup(
                lerpDup(fromx, tox, dx),
                lerpDup(fromy, toy, dx),
                dy,
            );
        }

        pub fn lerp(self: *Self, to: Self, d: eT) void {
            self.v = interp.lerp(
                Self.Vector4,
                .{ 1, 1, 1, 1 },
                self.v,
                to.v,
                Self.splat(d).v,
            );
        }
    };
}

pub const Vec4f = Vec4(f32);
pub const Vec4i = Vec4(i32);

pub fn Vec4fToVec4i(v: Vec4f) Vec4i {
    return .{ .v = .{
        @as(Vec4i.ElementType, @intFromFloat(v.v[0])),
        @as(Vec4i.ElementType, @intFromFloat(v.v[1])),
        @as(Vec4i.ElementType, @intFromFloat(v.v[2])),
        @as(Vec4i.ElementType, @intFromFloat(v.v[3])),
    } };
}

pub fn Vec4iToVec4f(v: Vec4i) Vec4f {
    return .{ .v = .{
        @as(Vec4f.ElementType, @floatFromInt(v.v[0])),
        @as(Vec4f.ElementType, @floatFromInt(v.v[1])),
        @as(Vec4f.ElementType, @floatFromInt(v.v[2])),
        @as(Vec4f.ElementType, @floatFromInt(v.v[3])),
    } };
}

// pub const Vec4f = struct {
//     v: Vector4f,

//     pub const Vector4f = @Vector(4, f32);

//     const Axis = enum(u8) {
//         X = 0,
//         Y = 1,
//         Z = 2,
//         W = 3,
//     };

//     pub inline fn initVector(v: Vector4f) Vec4f {
//         return .{ .v = v };
//     }

//     pub inline fn init(_x: f32, _y: f32, _z: f32, _w: f32) Vec4f {
//         return .{ .v = .{ _x, _y, _z, _w } };
//     }

//     pub inline fn splat(s: f32) Vec4f {
//         return .{ .v = .{ s, s, s, s } };
//     }

//     pub inline fn splat3(s: f32, _w: f32) Vec4f {
//         return .{ .v = .{ s, s, s, _w } };
//     }

//     pub fn zero() Vec4f {
//         return Vec4f.init(0, 0, 0, 0);
//     }

//     pub fn one() Vec4f {
//         return Vec4f.init(1, 1, 1, 1);
//     }

//     pub fn half() Vec4f {
//         return Vec4f.init(0.5, 0.5, 0.5, 0.5);
//     }

//     pub fn forward() Vec4f {
//         return Vec4f.init(0, 0, 1, 1);
//     }

//     pub fn up() Vec4f {
//         return Vec4f.init(0, 1, 0, 1);
//     }

//     pub fn right() Vec4f {
//         return Vec4f.init(1, 0, 0, 1);
//     }

//     // TODO: generate swizzle combinations?
//     pub inline fn x(self: Vec4f) f32 {
//         return self.v[@intFromEnum(Axis.X)];
//     }

//     pub inline fn y(self: Vec4f) f32 {
//         return self.v[@intFromEnum(Axis.Y)];
//     }

//     pub inline fn z(self: Vec4f) f32 {
//         return self.v[@intFromEnum(Axis.Z)];
//     }

//     pub inline fn w(self: Vec4f) f32 {
//         return self.v[@intFromEnum(Axis.W)];
//     }

//     pub inline fn setValue(self: *Vec4f, element: Axis, value: f32) void {
//         self.v[@intFromEnum(element)] = value;
//     }

//     pub inline fn setX(self: *Vec4f, value: f32) void {
//         self.setValue(Axis.X, value);
//     }

//     pub inline fn setY(self: *Vec4f, value: f32) void {
//         self.setValue(Axis.Y, value);
//     }

//     pub inline fn setZ(self: *Vec4f, value: f32) void {
//         self.setValue(Axis.Z, value);
//     }

//     pub inline fn setW(self: *Vec4f, value: f32) void {
//         self.setValue(Axis.W, value);
//     }

//     pub inline fn setXY(self: *Vec4f, sx: f32, sy: f32) void {
//         self.setValue(Axis.X, sx);
//         self.setValue(Axis.Y, sy);
//     }

//     pub inline fn setXYZ(self: *Vec4f, sx: f32, sy: f32, sz: f32) void {
//         self.setValue(Axis.X, sx);
//         self.setValue(Axis.Y, sy);
//         self.setValue(Axis.Z, sz);
//     }

//     pub inline fn fromXY(self: *Vec4f, sx: f32, sy: f32) void {
//         self.v = .{ sx, sy, 0.0, 1.0 };
//     }

//     pub inline fn fromXYZ(self: *Vec4f, sx: f32, sy: f32, sz: f32) void {
//         self.v = .{ sx, sy, sz, 0.0 };
//     }

//     pub inline fn set(self: *Vec4f, other: Vec4f) void {
//         self.v = other.v;
//     }

//     pub inline fn add(self: *Vec4f, other: Vec4f) void {
//         self.v += other.v;
//     }

//     pub inline fn addScalar(self: *Vec4f, scalar: f32) void {
//         self.v += @as(Vector4f, @splat(scalar));
//     }

//     pub inline fn addScalarDup(self: Vec4f, scalar: f32) Vec4f {
//         return Vec4f{ .v = self.v + @as(Vector4f, @splat(scalar)) };
//     }

//     pub inline fn addDup(self: Vec4f, other: Vec4f) Vec4f {
//         return Vec4f{ .v = self.v + other.v };
//     }

//     pub inline fn mul(self: *Vec4f, other: Vec4f) void {
//         self.v *= other.v;
//     }

//     pub inline fn mulDup(self: Vec4f, other: Vec4f) Vec4f {
//         return Vec4f{ .v = self.v * other.v };
//     }

//     pub inline fn scale(self: *Vec4f, scalar: f32) void {
//         self.v *= @as(Vector4f, @splat(scalar));
//     }

//     // pub inline fn scale3(self: *Vec4f, scalar: f32) void {
//     //     const vscalar = self.v;
//     //     vscalar[3] = 1.0;
//     //     self.v *= scalar;
//     // }

//     pub inline fn scaleDup(self: Vec4f, scalar: f32) Vec4f {
//         return Vec4f{ .v = self.v * @as(Vector4f, @splat(scalar)) };
//     }

//     pub inline fn scale3Dup(self: Vec4f, scalar: f32) Vec4f {
//         var v = self.v;
//         var vscalar = @as(Vector4f, @splat(scalar));
//         v[3] = 1.0;
//         vscalar[3] = 1.0;
//         return Vec4f{ .v = v * vscalar };
//     }

//     pub inline fn div(self: *Vec4f, scalar: f32) void {
//         self.v /= @as(Vector4f, @splat(scalar));
//     }

//     pub inline fn div3(self: *Vec4f, scalar: f32) void {
//         var v = self.v;
//         var vscalar = @as(Vector4f, @splat(scalar));
//         v[3] = 1.0;
//         vscalar[3] = 1.0;
//         self.v = v / vscalar;
//     }

//     pub inline fn divDup(self: Vec4f, scalar: f32) Vec4f {
//         return Vec4f{ .v = self.v / @as(Vector4f, @splat(scalar)) };
//     }

//     pub inline fn divVec(self: *Vec4f, other: Vec4f) void {
//         self.v /= other.v;
//     }

//     pub inline fn divVecDup(self: Vec4f, other: Vec4f) Vec4f {
//         return Vec4f{ .v = self.v / other.v };
//     }

//     pub inline fn dot3(self: Vec4f, other: Vec4f) f32 {
//         var vself = self.v;
//         vself[3] = 0;
//         var vother = other.v;
//         vother[3] = 0;
//         return @reduce(.Add, (vself * vother));
//         //return self.x() * other.x() + self.y() * other.y() + self.z() * other.z();
//     }

//     pub inline fn abs(self: Vec4f) Vec4f {
//         return .{ .v = @fabs(self.v) };
//     }

//     pub inline fn dot(self: Vec4f, other: Vec4f) f32 {
//         return @reduce(.Add, (self.v * other.v));
//     }

//     inline fn dotRaw(self: Vec4f, vother: Vector4f) f32 {
//         return @reduce(.Add, (self.v * vother));
//     }

//     pub inline fn cross3(self: Vec4f, other: Vec4f) Vec4f {
//         return Vec4f.init(self.y() * other.z() - self.z() * other.y(), self.z() * other.x() - self.x() * other.z(), self.x() * other.y() - self.y() * other.x(), 1);
//     }

//     pub inline fn sub(self: *Vec4f, other: Vec4f) void {
//         self.v -= other.v;
//     }

//     pub inline fn subDup(self: Vec4f, other: Vec4f) Vec4f {
//         return .{ .v = self.v - other.v };
//     }

//     pub inline fn subScalar(self: *Vec4f, scalar: f32) void {
//         self.v -= @as(Vector4f, @splat(scalar));
//     }

//     /// Returns vector length
//     pub inline fn length(self: Vec4f) f32 {
//         return math.sqrt(self.dot(self));
//     }

//     /// return length of vector squared. Avoids `math.sqrt`
//     pub inline fn lengthSqr(self: Vec4f) f32 {
//         return self.dot(self);
//     }

//     /// make `length` of vector 1.0 while maintaining direction
//     pub inline fn normalize(self: *Vec4f) void {
//         self.div(self.length());
//     }

//     /// constant version of normalize that returns a new `Vec4f` with length of 1.0
//     pub inline fn normalized(self: Vec4f) Vec4f {
//         return self.divDup(self.length());
//     }

//     pub inline fn length3(self: Vec4f) f32 {
//         return math.sqrt(self.length3Sqr());
//     }

//     /// return length of vector squared. Avoids `math.sqrt`
//     pub inline fn length3Sqr(self: Vec4f) f32 {
//         return self.x() * self.x() + self.y() * self.y() + self.z() * self.z();
//     }

//     /// make `length` of vector 1.0 while maintaining direction
//     pub inline fn normalize3(self: *Vec4f) void {
//         self.div3(self.length3());
//         // self.normaize();
//     }

//     /// constant version of normalize that returns a new `Vec4f` with length of 1.0
//     pub inline fn normalized3(self: Vec4f) Vec4f {
//         const len = self.length3();
//         return Vec4f.init(self.x() / len, self.y() / len, self.z() / len, self.w());
//         //return self.normalized();
//     }

//     pub inline fn clamp01(self: *Vec4f) void {
//         self.clamp(0.0, 1.0);
//     }

//     pub inline fn clamped01(self: Vec4f) Vec4f {
//         return .{ .v = self.clamped(0.0, 1.0) };
//     }

//     inline fn clamped(self: Vec4f, minVec4: f32, maxVec4: f32) Vector4f {
//         const maxV = @as(Vector4f, @splat(maxVec4));
//         const minV = @as(Vector4f, @splat(minVec4));

//         return @min(maxV, @max(minV, self.v));
//     }

//     pub inline fn maxScalar(self: Vec4f, scalar: f32) Vec4f {
//         const maxV = @as(Vector4f, @as(Vector4f, @splat(scalar)));
//         return .{ .v = @max(self.v, maxV) };
//     }

//     pub inline fn clamp(self: *Vec4f, minVec4: f32, maxVec4: f32) void {
//         self.v = self.clamped(minVec4, maxVec4);
//     }

//     inline fn clampedVec(self: Vec4f, minV: Vec4f, maxV: Vec4f) Vector4f {
//         return @min(maxV.v, @max(minV.v, self.v));
//     }

//     pub inline fn clampVec(self: *Vec4f, minVec4: Vec4f, maxVec4: Vec4f) void {
//         self.v = self.clampedVec(minVec4, maxVec4);
//     }

//     pub inline fn ceil(self: *Vec4f) void {
//         self.v = @ceil(self.v);
//     }

//     pub fn print(self: Vec4f) void {
//         std.debug.print(" [{} | {}]", .{ self.v, self.length() });
//     }

//     pub fn println(self: Vec4f) void {
//         self.print();
//         std.debug.print("\n", .{});
//     }

//     pub fn eq(self: Vec4f, other: Vec4f, _epsilon: f32) bool {
//         const vabs = @fabs(self.v - other.v);
//         return @reduce(.And, vabs < _epsilon);
//     }

//     pub inline fn min(lhs: Vec4f, rhs: Vec4f) Vec4f {
//         return .{ .v = @min(lhs.v, rhs.v) };
//     }

//     pub inline fn max(lhs: Vec4f, rhs: Vec4f) Vec4f {
//         return .{ .v = @max(lhs.v, rhs.v) };
//     }

//     pub inline fn maxElement(self: Vec4f) f32 {
//         return @reduce(.Max, self.v);
//     }

//     pub inline fn minElement(self: Vec4f) f32 {
//         return @reduce(.Min, self.v);
//     }

//     pub fn triArea(a: Vec4f, b: Vec4f, c: Vec4f) f32 {
//         return (c.x() - a.x()) * (b.y() - a.y()) -
//             (c.y() - a.y()) * (b.x() - a.x());
//     }

//     pub fn triBarycentericCoordsOld(v0: Vec4f, v1: Vec4f, v2: Vec4f, p: Vec4f) Vec4f {
//         return Vec4f.init(triArea(v1, v2, p), triArea(v2, v0, p), triArea(v0, v1, p), 0);
//     }

//     pub fn triBarycentericCoords(a: Vec4f, b: Vec4f, c: Vec4f, p: Vec4f) Vec4f {
//         const v0 = b.subDup(a);
//         const v1 = c.subDup(a);
//         const v2 = p.subDup(a);
//         const d00 = v0.dot(v0);
//         const d01 = v0.dot(v1);
//         const d11 = v1.dot(v1);
//         const d20 = v2.dot(v0);
//         const d21 = v2.dot(v1);
//         const invDenom = 1.0 / (d00 * d11 - d01 * d01);

//         const v = (d11 * d20 - d01 * d21) * invDenom;
//         const _w = (d00 * d21 - d01 * d20) * invDenom;
//         const u = 1.0 - v - _w;
//         return Vec4f.init(v, _w, u, 0.0);
//     }

//     pub fn triInterpStd(tri: Vec4f, v0: Vec4f, v1: Vec4f, v2: Vec4f, depth: f32, _w: f32) Vec4f {
//         return Vec4f.init((tri.x() * v0.x() + tri.y() * v1.x() + tri.z() * v2.x()) / depth, (tri.x() * v0.y() + tri.y() * v1.y() + tri.z() * v2.y()) / depth, (tri.x() * v0.z() + tri.y() * v1.z() + tri.z() * v2.z()) / depth, _w);
//     }

//     pub fn triInterp(tri: Vec4f, v0: Vec4f, v1: Vec4f, v2: Vec4f, depth: f32, _w: f32) Vec4f {
//         const vx = Vector4f{ v0.x(), v1.x(), v2.x(), 0.0 };
//         const vy = Vector4f{ v0.y(), v1.y(), v2.y(), 0.0 };
//         const vz = Vector4f{ v0.z(), v1.z(), v2.z(), 0.0 };

//         return Vec4f.init(tri.dotRaw(vx) / depth, tri.dotRaw(vy) / depth, tri.dotRaw(vz) / depth, _w);
//     }

//     pub fn triInterpArray(tri: Vec4f, v: [3]Vec4f, depth: f32, _w: f32) Vec4f {
//         return triInterp(tri, v[0], v[1], v[2], depth, _w);
//     }

//     pub fn lerpDup(from: Vec4f, to: Vec4f, d: f32) Vec4f {
//         return Vec4f.init(interp.lerp(f32, 1.0, from.x(), to.x(), d), interp.lerp(f32, 1.0, from.y(), to.y(), d), interp.lerp(f32, 1.0, from.z(), to.z(), d), interp.lerp(f32, 1.0, from.w(), to.w(), d));
//     }

//     pub fn blerpDup(fromx: Vec4f, tox: Vec4f, fromy: Vec4f, toy: Vec4f, dx: f32, dy: f32) Vec4f {
//         return lerpDup(lerpDup(fromx, tox, dx), lerpDup(fromy, toy, dx), dy);
//     }

//     pub fn lerp(self: *Vec4f, from: Vec4f, to: Vec4f, d: f32) void {
//         self.x = interp.lerp(1.0, from.x(), to.x(), d);
//         self.y = interp.lerp(1.0, from.y(), to.y(), d);
//         self.z = interp.lerp(1.0, from.z(), to.z(), d);
//         self.w = interp.lerp(1.0, from.w(), to.w(), d);
//     }

//     // pub fn eq3(self: Vec4f, other: Vec4f, _epsilon: f32) bool {
//     //     return (math.fabs(self.x - other.x) < _epsilon) and
//     //         (math.fabs(self.y - other.y) < _epsilon) and
//     //         (math.fabs(self.z - other.z) < _epsilon);
//     // }

//     // pub fn triArea(a: Vec4f, b: Vec4f, c: Vec4f) f32 {
//     //     return (c.x() - a.x()) *
//     //            (b.y() - a.y()) -
//     //            (c.y() - a.y()) *
//     //            (b.x() - a.x());
//     // }

//     // pub fn triCoords(v0: Vec4f, v1: Vec4f, v2: Vec4f, p: Vec4f) Vec4f {
//     //     return Vec4f.init(
//     //         triArea(v1, v2, p),
//     //         triArea(v2, v0, p),
//     //         triArea(v0, v1, p),
//     //         0);
//     // }

//     // pub fn triInterp(tri: Vec4f, v0: Vec4f, v1: Vec4f, v2: Vec4f, depth: f32, _w: f32) Vec4f {
//     //     var t = Vec4f{.v=tri.v * v0.v};
//     //     t.v += (tri.v * v1.v);
//     //     t.v += (tri.v * v2.v);
//     //     t.v /= @splat(depth);
//     //     t.v[3] = _w;

//     //     return t;
//     // }

//     // pub fn triInterpArray(tri: Vec4f, v: [3]Vec4f, depth: f32, _w: f32) Vec4f {
//     //     return triInterp(tri, v[0], v[1], v[2], depth, _w);
//     // }

//     // pub fn lerpDup(from: Vec4f, to: Vec4f, d: f32) Vec4f {
//     //     return .{.v=interp.lerp(Vector4f, @Vector(4,f32){1.0,1.0,1.0,1.0}, from.v, to.v, @splat( d))};
//     // }

//     // pub fn blerpDup(fromx: Vec4f, tox: Vec4f, fromy: Vec4f, toy: Vec4f, dx: f32, dy: f32) Vec4f {
//     //     return lerpDup(lerpDup(fromx, tox, dx), lerpDup(fromy, toy, dx), dy);
//     // }

//     // pub fn lerp(self: *Vec4f, from: Vec4f, to: Vec4f, d: f32) void {
//     //     self.v = interp.lerp(Vector4f, @Vector(4,f32){1.0,1.0,1.0,1.0}, from.v, to.v, d);
//     // }
// };

const assert = @import("std").debug.assert;

const epsilon = 0.00001;

pub fn assert_f32_equal(actual: f32, expected: f32) void {
    const diff = math.fabs(actual - expected);
    // std.debug.print("{} - {} = {}, e={}, d < e = {}\n", .{actual, expected, diff, epsilon, diff < epsilon });
    assert(diff < epsilon);
}

test "Vec4f.alignment" {
    assert(@sizeOf(Vec4f) == (@sizeOf(f32) * 4));
    assert(@alignOf(Vec4f) == @alignOf(@Vector(4, f32)));
}

test "Vec4f.add" {
    var lhs = Vec4f.init(1.0, 2.0, 3.0, 1.0);
    const rhs = Vec4f.init(2.0, 3.0, 4.0, 1.0);

    // std.debug.print("\n", .{});
    // lhs.println();
    // rhs.println();
    lhs.add(rhs);

    // std.debug.print("after add:\n", .{});
    // lhs.println();
    //rhs.println();

    assert_f32_equal(lhs.x(), 3.0);
    assert_f32_equal(lhs.y(), 5.0);
    assert_f32_equal(lhs.z(), 7.0);
    assert_f32_equal(lhs.w(), 2.0);
}

test "Vec4f.set" {
    var lhs = Vec4f.init(1.0, 2.0, 3.0, 1.0);
    const rhs = Vec4f.init(2.0, 3.0, 4.0, 1.0);
    lhs.set(rhs);

    assert_f32_equal(lhs.x(), rhs.x());
    assert_f32_equal(lhs.y(), rhs.y());
    assert_f32_equal(lhs.z(), rhs.z());
    assert_f32_equal(lhs.w(), rhs.w());
}

test "Vec4f.normalize" {
    var lhs = Vec4f.init(1.0, 2.0, 3.0, 1.0);
    lhs.normalize();
    assert_f32_equal(lhs.length(), 1.0);
}

test "Vec4f.normalize3" {
    var lhs = Vec4f.init(1.0, 2.0, 3.0, 1.0);
    lhs.normalize3();
    std.debug.print("length = {}, {}", .{ lhs.length3(), lhs.v });
    assert_f32_equal(lhs.length3(), 1.0);
    assert_f32_equal(lhs.w(), 1.0);
}

test "Vec4f.normalized" {
    const lhs = Vec4f.init(1.0, 2.0, 3.0, 1.0);
    const lhsLen = lhs.length();
    const normal = lhs.normalized();
    assert_f32_equal(normal.length(), 1.0);
    assert_f32_equal(lhs.length(), lhsLen);
}

test "Vec4f.lengthSqr" {
    const lhs = Vec4f.init(1.0, 2.0, 3.0, 1.0);
    const len = lhs.length();
    const sqr = lhs.lengthSqr();
    assert_f32_equal(sqr, len * len);
}
