const std = @import("std");
const warn = std.debug.print;
const fmt = std.fmt;
const assert = @import("std").debug.assert;
const math = std.math;
const io = std.io;
const Allocator = std.mem.Allocator;

const Vec4f = @import("../core/vector.zig").Vec4f;
const Mat44f = @import("../core/matrix.zig").Mat44f;
const Profile = @import("../core/profiler.zig").Profile;

const Mesh = @import("../../render/mesh.zig").Mesh;

pub const MeshObjLoader = @import("render/obj_mesh_loader.zig");
pub const TgaTexLoader = @import("render/tga_texture_loader.zig");
pub const R32TexLoader = @import("render/r32_texture_loader.zig");

pub const Importer = struct {
    pub fn importObjMesh(alloc: *Allocator, objfile: []const u8, assetfile: []const u8) !Mesh {}
};
