// pub const tracy = @import("tracy/tracy-cImport.zig");
pub const tracy = @import("tracy/tracy-zig.zig");
pub const trace = tracy.trace;
pub const traceNamed = tracy.traceNamed;
// pub const trace = @import("tracy/tracy-zig.zig").trace;
// pub const trace = @import("tracy/tracy-zig-srcloc.zig").trace;
