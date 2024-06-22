pub const Header = packed struct {
    marker: [4]u8,
    version: u32,

    index: Descriptor,
};

pub const Index = packed struct {
    start: u64,
    count: u32,
};

pub const Descriptor = packed struct {
    type: u32,
    offset: u64,
    size: u64,
};


pub fn write() void {

}