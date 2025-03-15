// Header
//  Descriptor
//  Data
//  ...
const std = @import("std");
const UUID = @import("../uuid.zig").UUID;

pub const Header = extern struct {
    marker: [8]u8 = [_]u8{ 'A', 'S', 'S', 'E', 'T', 0, 0, 0 },
    version: u32,
};

pub const Descriptor = extern struct {
    type: u32,
    offset: u64,
    size: u64,
    uuid: u128,
    hash: u128,
};

test "Write" {
    const header = Header{
        .version = 1,
    };

    try testing.expectEqual(1, header.version);
}

const ArrayList = std.ArrayList;
const test_allocator = std.testing.allocator;
const testing = std.testing;

test "io writer usage" {
    var list = ArrayList(u8).init(test_allocator);
    defer list.deinit();

    const header = Header{
        .version = 1,
    };

    const writer = list.writer();
    try writer.writeStruct(header);
    const bytes_written = list.items.len;

    try testing.expect(bytes_written == @sizeOf(Header));

    var readStream = std.io.fixedBufferStream(list.items);
    var reader = readStream.reader();

    for (list.items) |d| {
        std.debug.print("{X} ", .{d});
    }

    std.debug.print("\n{}\n", .{header});

    const readHeader = try reader.readStruct(Header);
    try testing.expectEqual(header, readHeader);
    std.debug.print("\n{}\n", .{readHeader});
}
