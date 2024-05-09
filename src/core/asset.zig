pub const Header = packed struct {
    marker: [4]u8,
};

pub const Descriptor = packed struct {
    type: u32,
    offset: u64,
    size: u64,
};
