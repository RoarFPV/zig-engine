const std = @import("std");
const math = std.math;
const vector = @import("vector.zig");
const Vec4f = vector.Vec4f;
const assert_f32_equal = vector.assert_f32_equal;

pub const PI = 3.1415926535897932384626433832795;
pub const DEG2RAD = PI / 180.0;

pub const Mat44f = struct {
    mm: [4]Vec4f,

    pub inline fn init(self: *Mat44f, values: [4]Vec4f) void {
        self.mm = values;
    }

    pub inline fn copy(self: *Mat44f, other: Mat44f) void {
        self.mm[0] = other.mm[0];
        self.mm[1] = other.mm[1];
        self.mm[2] = other.mm[2];
        self.mm[3] = other.mm[3];
    }

    pub inline fn identity() Mat44f {
        return Mat44f{
            .mm = [4]Vec4f{
                Vec4f.init(1, 0, 0, 0),
                Vec4f.init(0, 1, 0, 0),
                Vec4f.init(0, 0, 1, 0),
                Vec4f.init(0, 0, 0, 1),
            },
        };
    }

    pub inline fn rotX(radAngle: f32) Mat44f {
        const cosA = math.cos(radAngle);
        const sinA = math.sin(radAngle);
        return Mat44f{
            .mm = [4]Vec4f{
                Vec4f.init(1, 0, 0, 0),
                Vec4f.init(0, cosA, -sinA, 0),
                Vec4f.init(0, sinA, cosA, 0),
                Vec4f.init(0, 0, 0, 1),
            },
        };
    }

    pub inline fn rotY(radAngle: f32) Mat44f {
        const cosA = math.cos(radAngle);
        const sinA = math.sin(radAngle);
        return Mat44f{
            .mm = [4]Vec4f{
                Vec4f.init(cosA, 0, sinA, 0),
                Vec4f.init(0, 1, 0, 0),
                Vec4f.init(-sinA, 0, cosA, 0),
                Vec4f.init(0, 0, 0, 1),
            },
        };
    }

    pub inline fn rotZ(radAngle: f32) Mat44f {
        const cosA = math.cos(radAngle);
        const sinA = math.sin(radAngle);
        return Mat44f{
            .mm = [4]Vec4f{
                Vec4f.init(cosA, -sinA, 0, 0),
                Vec4f.init(sinA, cosA, 0, 0),
                Vec4f.init(0, 0, 1, 0),
                Vec4f.init(0, 0, 0, 1),
            },
        };
    }

    pub inline fn mul_vec4(self: Mat44f, vec: Vec4f) Vec4f {
        return Vec4f.init(
            self.mm[0].dot(vec),
            self.mm[1].dot(vec),
            self.mm[2].dot(vec),
            self.mm[3].dot(vec),
        );
    }

    pub inline fn col(self: Mat44f, index: u2) Vec4f {
        switch (index) {
            0 => return Vec4f.init(self.mm[0].x(), self.mm[1].x(), self.mm[2].x(), self.mm[3].x()),
            1 => return Vec4f.init(self.mm[0].y(), self.mm[1].y(), self.mm[2].y(), self.mm[3].y()),
            2 => return Vec4f.init(self.mm[0].z(), self.mm[1].z(), self.mm[2].z(), self.mm[3].z()),
            3 => return Vec4f.init(self.mm[0].w(), self.mm[1].w(), self.mm[2].w(), self.mm[3].w()),
            //else => @compileError("Invalid mat44 column index {}",.{index}),
        }
        unreachable;
    }

    pub fn mul(self: *Mat44f, other: Mat44f) void {
        self.mm = [4]Vec4f{
            Vec4f.init(self.mm[0].dot(other.col(0)), self.mm[0].dot(other.col(1)), self.mm[0].dot(other.col(2)), self.mm[0].dot(other.col(3))),
            Vec4f.init(self.mm[1].dot(other.col(0)), self.mm[1].dot(other.col(1)), self.mm[1].dot(other.col(2)), self.mm[1].dot(other.col(3))),
            Vec4f.init(self.mm[2].dot(other.col(0)), self.mm[2].dot(other.col(1)), self.mm[2].dot(other.col(2)), self.mm[2].dot(other.col(3))),
            Vec4f.init(self.mm[3].dot(other.col(0)), self.mm[3].dot(other.col(1)), self.mm[3].dot(other.col(2)), self.mm[3].dot(other.col(3))),
        };
    }

    pub inline fn mulDup(self: *const Mat44f, other: Mat44f) Mat44f {
        return Mat44f{ .mm = [4]Vec4f{
            Vec4f.init(self.mm[0].dot(other.col(0)), self.mm[0].dot(other.col(1)), self.mm[0].dot(other.col(2)), self.mm[0].dot(other.col(3))),
            Vec4f.init(self.mm[1].dot(other.col(0)), self.mm[1].dot(other.col(1)), self.mm[1].dot(other.col(2)), self.mm[1].dot(other.col(3))),
            Vec4f.init(self.mm[2].dot(other.col(0)), self.mm[2].dot(other.col(1)), self.mm[2].dot(other.col(2)), self.mm[2].dot(other.col(3))),
            Vec4f.init(self.mm[3].dot(other.col(0)), self.mm[3].dot(other.col(1)), self.mm[3].dot(other.col(2)), self.mm[3].dot(other.col(3))),
        } };
    }

    pub inline fn mul33(self: *Mat44f, other: Mat44f) void {
        self.mm = [4]Vec4f{
            Vec4f.init(self.mm[0].dot(other.col(0)), self.mm[0].dot(other.col(1)), self.mm[0].dot(other.col(2)), 0),
            Vec4f.init(self.mm[1].dot(other.col(0)), self.mm[1].dot(other.col(1)), self.mm[1].dot(other.col(2)), 0),
            Vec4f.init(self.mm[2].dot(other.col(0)), self.mm[2].dot(other.col(1)), self.mm[2].dot(other.col(2)), 0),
            self.mm[3],
        };
    }

    pub fn mul33_vec4(self: *const Mat44f, v: Vec4f) Vec4f {
        return Vec4f.init(
            self.mm[0].dot3(v),
            self.mm[1].dot3(v),
            self.mm[2].dot3(v),
            v.w(),
        );
    }

    pub fn mul33_divW_vec4(self: *const Mat44f, v: Vec4f) Vec4f {
        const w = self.mm[3].dot3(v);
        const invW = 1.0 / w;
        return Vec4f.init(
            self.mm[0].dot3(v) * invW,
            self.mm[1].dot3(v) * invW,
            self.mm[2].dot3(v) * invW,
            w,
        );
    }

    pub fn forward(self: Mat44f) Vec4f {
        return self.mul33_vec4(Vec4f.forward()).normalized();
    }

    pub inline fn up(self: Mat44f) Vec4f {
        return self.mul33_vec4(Vec4f.up()).normalized();
    }

    pub inline fn right(self: Mat44f) Vec4f {
        return self.mul33_vec4(Vec4f.right()).normalized();
    }

    pub fn createPerspective(fovy: f32, aspect: f32, nearZ: f32, farZ: f32) Mat44f {
        const tangent = math.tan(fovy / 2 * DEG2RAD); // tangent of half fovY
        const height = nearZ * tangent; // half height of near plane
        const width = height * aspect; // half width of near plane

        const left = -width;
        const r = width;
        const bottom = -height;
        const top = height;

        const deltaX = r - left;
        const deltaY = top - bottom;
        const deltaZ = farZ - nearZ;

        // if ((nearZ <= 0.0f) || (farZ <= 0.0f) || (deltaX <= 0.0f) || (deltaY <= 0.0f) || (deltaZ <= 0.0f))
        //   return;

        return Mat44f{
            .mm = [4]Vec4f{
                Vec4f.init((2.0 * nearZ) / deltaX, 0.0, (r + left) / deltaX, 0.0),
                Vec4f.init(0.0, (2.0 * nearZ) / deltaY, (top + bottom) / deltaY, 0.0),
                Vec4f.init(0.0, 0.0, -(farZ + nearZ) / deltaZ, (-2.0 * nearZ * farZ) / deltaZ),
                Vec4f.init(0.0, 0.0, -1.0, 0.0),
            },
        };
    }

    pub fn createPerspectiveSimple(fovy: f32, aspect: f32, nearZ: f32, farZ: f32) Mat44f {
        const tangent = math.tan(fovy / 2 * DEG2RAD);
        const deltaZ = nearZ - farZ;
        return Mat44f{
            .mm = [4]Vec4f{
                Vec4f.init(1.0 / (tangent * aspect), 0.0, 0.0, 0.0),
                Vec4f.init(0.0, 1.0 / tangent, 0.0, 0.0),
                Vec4f.init(0.0, 0.0, (-farZ - nearZ) / deltaZ, 2.0 * nearZ * farZ / deltaZ),
                // Vec4f.init(0.0, 0.0, 1.0, 0.0),
                Vec4f.init(0.0, 0.0, 1.0, 0.0),
            },
        };
    }

    pub fn translate(self: *Mat44f, vec: Vec4f) void {
        const tmp = Mat44f{
            .mm = [4]Vec4f{
                Vec4f.init(1, 0, 0, vec.x()),
                Vec4f.init(0, 1, 0, vec.y()),
                Vec4f.init(0, 0, 1, vec.z()),
                Vec4f.init(0, 0, 0, 1),
            },
        };

        self.mul(tmp);
    }

    pub fn position(self: Mat44f) Vec4f {
        return .{ .v = .{ self.mm[0].w(), self.mm[1].w(), self.mm[2].w(), 0.0 } };
    }

    pub fn print(self: Mat44f) void {
        std.debug.print("\n[\n", .{});
        self.mm[0].println();
        self.mm[1].println();
        self.mm[2].println();
        self.mm[3].println();
        std.debug.print("]\n", .{});
    }
};

test "Mat44f.mul" {
    var a = Mat44f.identity();
    const b = Mat44f.identity();

    std.debug.print("\na:", .{});
    a.print();

    std.debug.print("\nb:", .{});
    b.print();

    a.mul(b);

    std.debug.print("\na result:", .{});
    a.print();
}

test "Mat44f.translate" {
    const lhs = Vec4f.init(2.0, 2.0, 2.0, 1.0);
    var model = Mat44f.identity();

    model.print();
    model.translate(lhs);
    model.print();

    const origin = model.mul_vec4(Vec4f.init(0, 0, 0, 1));
    lhs.print();
    origin.print();
    //const mat = Mat44f.init([0.0]f32 ** (4*4));

    assert_f32_equal(lhs.x(), origin.x());
    assert_f32_equal(lhs.y(), origin.y());
    std.debug.print("\n{} ? {}\n", .{ lhs.z(), origin.z() });
    assert_f32_equal(lhs.z(), origin.z());
    assert_f32_equal(lhs.w(), origin.w());
}
