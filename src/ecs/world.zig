const std = @import("std");

pub const Entity = struct {
    const Self = @This();

    id: usize = 0,
};

pub const World = struct {
    const Self = @This();

    const EntityMap = std.AutoArrayHashMap(usize, Entity);

    pub fn init() Self {
        return .{};
    }
};
