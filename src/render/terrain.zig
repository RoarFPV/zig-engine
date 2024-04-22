const std = @import("std");
const warn = std.debug.print;
const fmt = std.fmt;
const assert = @import("std").debug.assert;
const math = std.math;
const io = std.io;

const vector = @import("../core/vector.zig");
const Vec4f = vector.Vec4f;
const Mat44f = @import("../core/matrix.zig").Mat44f;
const Profile = @import("../core/profiler.zig").Profile;
const Bounds = @import("../core/bounds.zig").Bounds;

const Texture = @import("texture.zig").Texture;

pub const Terrain = struct {
    vertexBuffer: []Vec4f,
    // vertexNormalBuffer: []Vec4f,
    // textureCoordBuffer: []Vec4f,
    // colorBuffer: []Vec4f,
    bounds: Bounds,
    size: vector.Vec4i,

    pub fn init(allocator: *std.mem.Allocator, heightmap: Texture, width: usize, height: usize) !Terrain {
        assert(heightmap.widthBytes == width);
        assert(heightmap.heightBytes == height);
        const total = width * height;
        var bounds = Bounds.init(Vec4f.zero(), Vec4f.zero());
        var verts = try allocator.alloc(Vec4f, total);

        for (0..total) |index| {
            const vertHeight = heightmap.sampleIndex(index);
            const y = index / width;
            const x = index - y * width;

            const v = Vec4f.init(@floatFromInt(x), vertHeight.x(), @floatFromInt(y), 1);
            bounds.add(v);
            verts[index] = v;
        }

        return .{
            .vertexBuffer = verts,
            .bounds = bounds,
            .size = vector.Vec4i.init(@intCast(width), 1, @intCast(height), 0),
        };
    }
};

// pub fn recalculateNormals(self: *Mesh) void {
//     const ids = self.indexBuffer.len;
//     const numTris = ids / 3;
//     var normals = self.triNormals[0..];

//     var tri: u16 = 0;
//     while (tri < ids) {
//         const vi0 = self.indexBuffer[tri + 0];
//         const vi1 = self.indexBuffer[tri + 1];
//         const vi2 = self.indexBuffer[tri + 2];

//         const v0 = self.vertexBuffer[vi0];
//         const v1 = self.vertexBuffer[vi1];
//         const v2 = self.vertexBuffer[vi2];

//         var e0 = v1;
//         var e1 = v2;

//         e0.sub(v0);
//         e1.sub(v0);

//         const rnormal = e0.cross3(e1);
//         const normal = rnormal.normalized3();

//         normals[tri] = normal;
//         tri += 3;
//     }

//     var vn: u16 = 0;
//     while (vn < self.vertexNormalBuffer.len) {
//         defer vn += 1;
//         _ = self.vertexNormalBuffer[vn];

//         tri = 0;
//         while (tri < ids) {
//             tri += 1;

//             const v0 = self.vertexBuffer[self.indexBuffer[tri]];
//             _ = v0;
//         }
//     }
// }

// pub fn recalculateNormals_old(self: *Mesh) void {
//     const ids = self.indexBuffer.len;
//     // const numTris = ids / 3;

//     var tri: u16 = 0;
//     while (tri < ids) {
//         const vi0 = self.indexBuffer[tri + 0];
//         const vi1 = self.indexBuffer[tri + 1];
//         const vi2 = self.indexBuffer[tri + 2];

//         const v0 = self.vertexBuffer[vi0];
//         const v1 = self.vertexBuffer[vi1];
//         const v2 = self.vertexBuffer[vi2];

//         var e0 = v1;
//         var e1 = v2;

//         e0.sub(v0);
//         e1.sub(v0);

//         const rnormal = e0.cross3(e1);
//         const normal = rnormal.normalized3();

//         //warn("[{}:{}] : ", .{tri, vi0});
//         self.vertexNormalBuffer[vi0].add(normal);
//         self.vertexNormalBuffer[vi0].normalize3();

//         self.vertexNormalBuffer[vi1].add(normal);
//         self.vertexNormalBuffer[vi1].normalize3();

//         self.vertexNormalBuffer[vi2].add(normal);
//         self.vertexNormalBuffer[vi2].normalize3();

//         tri += 3;
//     }

//     tri = 0;
//     while (tri < self.vertexNormalBuffer.len) {
//         // warn("\n[{}] : ", .{tri});
//         self.vertexNormalBuffer[tri].w = 1;
//         self.vertexNormalBuffer[tri].normalize3();
//         tri += 1;
//     }
// }
