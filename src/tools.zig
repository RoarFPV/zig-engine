// platform imports
const std = @import("std");
const fmt = std.fmt;
const warn = std.debug.print;
const assert = std.debug.assert;
const Timer = std.time.Timer;

pub const MeshObjLoader = @import("tools/render/obj_mesh_loader.zig");
pub const TgaTexLoader = @import("tools/render/tga_texture_loader.zig");
pub const R32TexLoader = @import("tools/render/r32_texture_loader.zig");

pub var stdout = std.io.getStdOut();
