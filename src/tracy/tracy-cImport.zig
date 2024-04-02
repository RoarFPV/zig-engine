//! Tracy bindings from Zig compiler
//
// The MIT License (Expat)
//
// Copyright (c) 2015-2022, Zig contributors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

const std = @import("std");
const builtin = @import("builtin");
const build_options = @import("build_options");

const c = @cImport({
    @cDefine("TRACY_ENABLE", "1");
    //@cDefine("TRACY_FIBERS", "1");
    @cInclude("tracy/TracyC.h");
});

pub const enable = if (builtin.is_test) false else build_options.enable_tracy;
pub const enable_allocation = enable and build_options.enable_tracy_allocation;
pub const enable_callstack = enable and build_options.enable_tracy_callstack;

// TODO: make this configurable
const callstack_depth = 10;

const CtxWrapper = struct {
    ctx: ?c.___tracy_c_zone_context,

    pub fn init(_ctx: ?c.___tracy_c_zone_context) Ctx {
        return .{ .ctx = _ctx };
    }

    pub inline fn end(self: @This()) void {
        c.___tracy_emit_zone_end(self.ctx.?);
    }

    pub inline fn addText(self: @This(), text: []const u8) void {
        c.___tracy_emit_zone_text(self.ctx.?, text.ptr, text.len);
    }

    pub inline fn setName(self: @This(), name: []const u8) void {
        c.___tracy_emit_zone_name(self.ctx.?, name.ptr, name.len);
    }

    pub inline fn setColor(self: @This(), color: u32) void {
        c.___tracy_emit_zone_color(self.ctx.?, color);
    }

    pub inline fn setValue(self: @This(), value: u64) void {
        c.___tracy_emit_zone_value(self.ctx.?, value);
    }
};

pub const Ctx = if (enable) CtxWrapper else struct {
    pub inline fn end(self: @This()) void {
        _ = self;
    }

    pub inline fn addText(self: @This(), text: []const u8) void {
        _ = self;
        _ = text;
    }

    pub inline fn setName(self: @This(), name: []const u8) void {
        _ = self;
        _ = name;
    }

    pub inline fn setColor(self: @This(), color: u32) void {
        _ = self;
        _ = color;
    }

    pub inline fn setValue(self: @This(), value: u64) void {
        _ = self;
        _ = value;
    }
};

pub inline fn trace(comptime src: std.builtin.SourceLocation) Ctx {
    if (!enable) return .{};

    // TODO: the below `.line = 1,` should be `.line = src.line`, this is blocked by
    //       https://github.com/ziglang/zig/issues/13315

    if (enable_callstack) {
        return Ctx.init(c.___tracy_emit_zone_begin_callstack(c.___tracy_alloc_srcloc(
            src.line,
            src.file.ptr,
            src.file.len,
            src.fn_name.ptr,
            src.fn_name.len,
        ), callstack_depth, 1));
    } else {
        return Ctx.init(c.___tracy_emit_zone_begin(c.___tracy_alloc_srcloc(
            src.line,
            src.file.ptr,
            src.file.len,
            src.fn_name.ptr,
            src.fn_name.len,
        ), 1));
    }
}

pub inline fn traceNamed(comptime src: std.builtin.SourceLocation, comptime name: [:0]const u8) Ctx {
    if (!enable) return .{};

    // TODO: the below `.line = 1,` should be `.line = src.line`, this is blocked by
    //       https://github.com/ziglang/zig/issues/13315

    if (enable_callstack) {
        return Ctx.init(c.___tracy_emit_zone_begin_callstack(c.___tracy_alloc_srcloc_name(
            src.line,
            src.file.ptr,
            src.file.len,
            src.fn_name.ptr,
            src.fn_name.len,
            name.ptr,
            name.len,
        ), callstack_depth, 1));
    } else {
        return Ctx.init(c.___tracy_emit_zone_begin(c.___tracy_alloc_srcloc_name(
            src.line,
            src.file.ptr,
            src.file.len,
            src.fn_name.ptr,
            src.fn_name.len,
            name.ptr,
            name.len,
        ), 1));
    }
}

pub fn tracyAllocator(allocator: std.mem.Allocator) TracyAllocator(null) {
    return TracyAllocator(null).init(allocator);
}

pub fn TracyAllocator(comptime name: ?[:0]const u8) type {
    return struct {
        parent_allocator: std.mem.Allocator,

        const Self = @This();

        pub fn init(parent_allocator: std.mem.Allocator) Self {
            return .{
                .parent_allocator = parent_allocator,
            };
        }

        pub fn allocator(self: *Self) std.mem.Allocator {
            return .{
                .ptr = self,
                .vtable = &.{
                    .alloc = allocFn,
                    .resize = resizeFn,
                    .free = freeFn,
                },
            };
        }

        fn allocFn(ptr: *anyopaque, len: usize, ptr_align: u8, ret_addr: usize) ?[*]u8 {
            const self = @as(*Self, @ptrCast(@alignCast(@alignOf(Self), ptr)));
            const result = self.parent_allocator.rawAlloc(len, ptr_align, ret_addr);
            if (result) |data| {
                if (len != 0) {
                    if (name) |n| {
                        allocNamed(data, len, n);
                    } else {
                        alloc(data, len);
                    }
                }
            } else {
                messageColor("allocation failed", 0xFF0000);
            }
            return result;
        }

        fn resizeFn(ptr: *anyopaque, buf: []u8, buf_align: u8, new_len: usize, ret_addr: usize) bool {
            const self = @as(*Self, @ptrCast(@alignCast(@alignOf(Self), ptr)));
            if (self.parent_allocator.rawResize(buf, buf_align, new_len, ret_addr)) {
                if (name) |n| {
                    freeNamed(buf.ptr, n);
                    allocNamed(buf.ptr, new_len, n);
                } else {
                    free(buf.ptr);
                    alloc(buf.ptr, new_len);
                }

                return true;
            }

            // during normal operation the compiler hits this case thousands of times due to this
            // emitting messages for it is both slow and causes clutter
            return false;
        }

        fn freeFn(ptr: *anyopaque, buf: []u8, buf_align: u8, ret_addr: usize) void {
            const self = @as(*Self, @ptrCast(@alignCast(@alignOf(Self), ptr)));
            self.parent_allocator.rawFree(buf, buf_align, ret_addr);
            // this condition is to handle free being called on an empty slice that was never even allocated
            // example case: `std.process.getSelfExeSharedLibPaths` can return `&[_][:0]u8{}`
            if (buf.len != 0) {
                if (name) |n| {
                    freeNamed(buf.ptr, n);
                } else {
                    free(buf.ptr);
                }
            }
        }
    };
}

// This function only accepts comptime known strings, see `messageCopy` for runtime strings
pub inline fn message(comptime msg: [:0]const u8) void {
    if (!enable) return;
    c.___tracy_emit_messageL(msg.ptr, if (enable_callstack) callstack_depth else 0);
}

// This function only accepts comptime known strings, see `messageColorCopy` for runtime strings
pub inline fn messageColor(comptime msg: [:0]const u8, color: u32) void {
    if (!enable) return;
    c.___tracy_emit_messageLC(msg.ptr, color, if (enable_callstack) callstack_depth else 0);
}

pub inline fn messageCopy(msg: []const u8) void {
    if (!enable) return;
    c.___tracy_emit_message(msg.ptr, msg.len, if (enable_callstack) callstack_depth else 0);
}

pub inline fn messageColorCopy(msg: [:0]const u8, color: u32) void {
    if (!enable) return;
    c.___tracy_emit_messageC(msg.ptr, msg.len, color, if (enable_callstack) callstack_depth else 0);
}

pub inline fn frameMark() void {
    if (!enable) return;
    c.___tracy_emit_frame_mark(null);
}

pub inline fn frameMarkNamed(comptime name: [:0]const u8) void {
    if (!enable) return;
    c.___tracy_emit_frame_mark(name.ptr);
}

pub inline fn namedFrame(comptime name: [:0]const u8) Frame(name) {
    frameMarkStart(name);
    return .{};
}

pub fn Frame(comptime name: [:0]const u8) type {
    return struct {
        pub fn end(_: @This()) void {
            frameMarkEnd(name);
        }
    };
}

inline fn frameMarkStart(comptime name: [:0]const u8) void {
    if (!enable) return;
    c.___tracy_emit_frame_mark_start(name.ptr);
}

inline fn frameMarkEnd(comptime name: [:0]const u8) void {
    if (!enable) return;
    c.___tracy_emit_frame_mark_end(name.ptr);
}

// extern fn c.___tracy_emit_frame_mark_start(name: [*:0]const u8) void;
// extern fn c.___tracy_emit_frame_mark_end(name: [*:0]const u8) void;

inline fn alloc(ptr: [*]u8, len: usize) void {
    if (!enable) return;

    if (enable_callstack) {
        c.___tracy_emit_memory_alloc_callstack(ptr, len, callstack_depth, 0);
    } else {
        c.___tracy_emit_memory_alloc(ptr, len, 0);
    }
}

inline fn allocNamed(ptr: [*]u8, len: usize, comptime name: [:0]const u8) void {
    if (!enable) return;

    if (enable_callstack) {
        c.___tracy_emit_memory_alloc_callstack_named(ptr, len, callstack_depth, 0, name.ptr);
    } else {
        c.___tracy_emit_memory_alloc_named(ptr, len, 0, name.ptr);
    }
}

inline fn free(ptr: [*]u8) void {
    if (!enable) return;

    if (enable_callstack) {
        c.___tracy_emit_memory_free_callstack(ptr, callstack_depth, 0);
    } else {
        c.___tracy_emit_memory_free(ptr, 0);
    }
}

inline fn freeNamed(ptr: [*]u8, comptime name: [:0]const u8) void {
    if (!enable) return;

    if (enable_callstack) {
        c.___tracy_emit_memory_free_callstack_named(ptr, callstack_depth, 0, name.ptr);
    } else {
        c.___tracy_emit_memory_free_named(ptr, 0, name.ptr);
    }
}

// extern fn c.___tracy_emit_zone_begin(
//     srcloc: *const c.___tracy_source_location_data,
//     active: c_int,
// ) c.___tracy_c_zone_context;
// extern fn c.___tracy_emit_zone_begin_callstack(
//     srcloc: *const c.___tracy_source_location_data,
//     depth: c_int,
//     active: c_int,
//) c.___tracy_c_zone_context;
// extern fn c.___tracy_emit_zone_text(ctx: c.___tracy_c_zone_context, txt: [*]const u8, size: usize) void;
// extern fn c.___tracy_emit_zone_name(ctx: c.___tracy_c_zone_context, txt: [*]const u8, size: usize) void;
// extern fn c.___tracy_emit_zone_color(ctx: c.___tracy_c_zone_context, color: u32) void;
// extern fn c.___tracy_emit_zone_value(ctx: c.___tracy_c_zone_context, value: u64) void;
// extern fn c.___tracy_emit_zone_end(ctx: c.___tracy_c_zone_context) void;
// extern fn c.___tracy_emit_memory_alloc(ptr: *const anyopaque, size: usize, secure: c_int) void;
// extern fn c.___tracy_emit_memory_alloc_callstack(ptr: *const anyopaque, size: usize, depth: c_int, secure: c_int) void;
// extern fn c.___tracy_emit_memory_free(ptr: *const anyopaque, secure: c_int) void;
// extern fn c.___tracy_emit_memory_free_callstack(ptr: *const anyopaque, depth: c_int, secure: c_int) void;
// extern fn c.___tracy_emit_memory_alloc_named(ptr: *const anyopaque, size: usize, secure: c_int, name: [*:0]const u8) void;
// extern fn c.___tracy_emit_memory_alloc_callstack_named(ptr: *const anyopaque, size: usize, depth: c_int, secure: c_int, name: [*:0]const u8) void;
// extern fn c.___tracy_emit_memory_free_named(ptr: *const anyopaque, secure: c_int, name: [*:0]const u8) void;
// extern fn c.___tracy_emit_memory_free_callstack_named(ptr: *const anyopaque, depth: c_int, secure: c_int, name: [*:0]const u8) void;
// extern fn c.___tracy_emit_message(txt: [*]const u8, size: usize, callstack: c_int) void;
// extern fn c.___tracy_emit_messageL(txt: [*:0]const u8, callstack: c_int) void;
// extern fn c.___tracy_emit_messageC(txt: [*]const u8, size: usize, color: u32, callstack: c_int) void;
// extern fn c.___tracy_emit_messageLC(txt: [*:0]const u8, color: u32, callstack: c_int) void;
// extern fn c.___tracy_emit_frame_mark(name: ?[*:0]const u8) void;

// extern fn c.___tracy_alloc_srcloc( line:u32, file: ?[*:0]const u8, fileSz:usize, function:[*:0]const u8, functionSz:usize ) *const c.___tracy_source_location_data;
// extern fn c.___tracy_alloc_srcloc_name( line:u32, file: ?[*:0]const u8, fileSz:usize, function:[*:0]const u8, functionSz:usize, name: ?[*:0]const u8, nameSz:usize ) *const c.___tracy_source_location_data;

// const c.___tracy_source_location_data = extern struct {
//     name: ?[*:0]const u8,
//     function: [*:0]const u8,
//     file: [*:0]const u8,
//     line: u32,
//     color: u32,
// };
