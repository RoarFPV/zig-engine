// Header
//  Descriptor
//  Data
//  ...
const std = @import("std");
const UUID = @import("../uuid.zig").UUID;
const serialized = @import("serialized.zig");
const Descriptor = serialized.Descriptor;

pub const Index = struct {
    const Self = @This();
    const PersistId = UUID.init().toU128();

    pub const DescriptorMap = std.AutoArrayHashMap(UUID, Descriptor);

    pub const Errors = error{
        InvalidIndexDescriptor,
    };

    descriptors: std.AutoArrayHashMap(UUID, Descriptor),
    uuid: UUID,

    pub fn init(allocator: std.mem.Allocator) Index {
        return Index{
            .uuid = UUID.init(),
            .descriptors = DescriptorMap.init(allocator),
        };
    }

    pub fn read(
        reader: std.io.Reader,
        allocator: std.mem.Allocator,
    ) !Index {
        const desc = reader.readStruct(Descriptor);

        if (desc.type != PersistId) {
            return Errors.InvalidIndexDescriptor;
        }

        var index = Index{
            .uuid = desc.uuid,
            .descriptors = DescriptorMap.init(allocator),
        };

        try index.readInto(desc, reader);

        return index;
    }

    pub fn readInto(self: Index, desc: Descriptor, reader: std.io.Reader) !void {
        const count = desc.size / @sizeOf(Descriptor);
        self.descriptors.ensureTotalCapacity(count);
        for (0..count) |_| {
            const assetDesc = try reader.readStruct(Descriptor);
            self.descriptors.putAssumeCapacity(
                assetDesc.uuid,
                assetDesc,
            );
        }
    }

    pub fn write(self: Index, offset: u64, writer: std.io.Writer) !void {
        const desc = Descriptor{
            .type = PersistId,
            .offset = offset,
            .size = self.descriptors.count() / @sizeOf(Descriptor),
            .uuid = self.uuid,
        };
        writer.writeStruct(desc);

        var iter = self.descriptors.iterator();
        while (iter.next()) |item| {
            writer.writeStruct(item.value_ptr.*);
        }
    }
};

// ======== Tests ============
const ArrayList = std.ArrayList;
const test_allocator = std.testing.allocator;
const testing = std.testing;

test "io writer usage" {
    // var list = ArrayList(u8).init(test_allocator);
    // defer list.deinit();

    // const writer = list.writer();
    // try writer.writeStruct(header);
    // const bytes_written = list.items.len;

    // try testing.expect(bytes_written == @sizeOf(Header));

    // var readStream = std.io.fixedBufferStream(list.items);
    // var reader = readStream.reader();

    // for (list.items) |d| {
    //     std.debug.print("{X} ", .{d});
    // }

    // std.debug.print("\n{}\n", .{header});

    // const readHeader = try reader.readStruct(Header);
    // try testing.expectEqual(header, readHeader);
    // std.debug.print("\n{}\n", .{readHeader});
}
