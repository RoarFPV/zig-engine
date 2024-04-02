const std = @import("std");
const Allocator = std.mem.Allocator;
const Vec4f = @import("../engine.zig").Vec4f;

pub const World = struct {
    const Self = @This();
    pub const Location = Vec4f;
    const NodeMap = std.AutoArrayHashMap(@Vector(4, i32), u32);

    allocator: Allocator,
    nodes: NodeMap,

    pub fn init(allocator: Allocator) Self {
        return .{
            .allocator = allocator,
            .nodes = NodeMap.init(allocator),
        };
    }

    pub fn deinit(self: *Self) void {
        self.nodes.deinit(self.allocator);
    }
};
