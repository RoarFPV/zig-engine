const std = @import("std");
const warn = std.debug.print;
const fmt = std.fmt;
const assert = @import("std").debug.assert;
const math = std.math;
const io = std.io;
const Allocator = std.mem.Allocator;
const Texture = @import("../../render/texture.zig").Texture;
const TextureFormat = @import("../../render/texture.zig").Format;
const Vec4f = @import("../../core/vector.zig").Vec4f;

pub const Error = error{InvalidImageFormat};

pub fn importR32File(allocator: *Allocator, file_path: []const u8) !Texture {
    const cwd = std.fs.cwd();

    const resolvedPath = try std.fs.path.resolve(allocator.*, &[_][]const u8{file_path});
    defer allocator.free(resolvedPath);

    std.debug.print("path: {s}", .{resolvedPath});

    var file = try cwd.openFile(resolvedPath, .{});
    defer file.close();

    var stream_source = io.StreamSource{ .file = file };
    var in = stream_source.reader();

    const sizeBytes = try file.getEndPos();
    const points = sizeBytes / 4;

    // assumes square
    const size = std.math.sqrt(points);

    var data = try std.ArrayList(Vec4f).initCapacity(allocator.*, points);

    for (0..points) |_| {
        const v = try in.readIntLittle(u32);

        data.appendAssumeCapacity(Vec4f.splat3(@bitCast(v), 1));
    }

    const t = try Texture.initVec(
        allocator,
        TextureFormat.F32,
        size,
        size,
        try data.toOwnedSlice(),
    );

    return t;
}
