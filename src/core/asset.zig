// Header
//  Descriptor
//  Data
//  ...
const std = @import("std");
const UUID = @import("uuid.zig").UUID;
const serialized = @import("asset/serialized.zig");
const Descriptor = serialized.Descriptor;
const Header = serialized.Header;
const Index = @import("asset/index.zig").Index;
// ======== Tests ============

const ArrayList = std.ArrayList;
const test_allocator = std.testing.allocator;
const testing = std.testing;

test "io writer usage" {
    var list = ArrayList(u8).init(test_allocator);
    defer list.deinit();

    const header = serialized.Header{
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
