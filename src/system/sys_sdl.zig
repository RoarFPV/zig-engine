const std = @import("std");
const warn = std.debug.print;
const fmt = std.fmt;

const common = @import("sys_common.zig");
const input = @import("sys_input.zig");
const dvui = @import("dvui");

const SDLBackend = dvui.backend;
comptime {
    std.debug.assert(@hasDecl(SDLBackend, "SDLBackend"));
}

const c = SDLBackend.c;

const assert = @import("std").debug.assert;

pub const trace = @import("../tracy.zig").trace;

const SDL_WINDOWPOS_UNDEFINED = @as(c_int, @bitCast(c.SDL_WINDOWPOS_UNDEFINED_MASK));

const SDL_INIT_EVERYTHING =
    c.SDL_INIT_TIMER |
    c.SDL_INIT_AUDIO |
    c.SDL_INIT_VIDEO |
    c.SDL_INIT_EVENTS |
    c.SDL_INIT_JOYSTICK |
    c.SDL_INIT_HAPTIC |
    c.SDL_INIT_GAMECONTROLLER;

var t0: u32 = 0;
//var t1:u32 = 0;

pub const Config = struct {
    windowWidth: u16,
    windowHeight: u16,

    renderWidth: u16,
    renderHeight: u16,

    maxFps: u16,
    fullscreen: bool,

    pub inline fn targetDt(self: Config) u16 {
        _ = self;
        return 1000 / config.maxFps;
    }
};

var renderTexture: ?*c.SDL_Texture = null;
var window: ?*c.SDL_Window = null;
var renderer: ?*c.SDL_Renderer = null;

var uiWindow: dvui.Window = undefined;
var uiSDLBackend: SDLBackend = undefined;

var gpa_instance = std.heap.GeneralPurposeAllocator(.{}){};
const gpa = gpa_instance.allocator();

var config = Config{
    .windowWidth = 800,
    .windowHeight = 600,
    .renderWidth = 800,
    .renderHeight = 600,
    .maxFps = 60,
    .fullscreen = false,
};

pub fn init(cfg: Config) !void {
    config = cfg;

    if (c.SDL_Init(c.SDL_INIT_TIMER | c.SDL_INIT_AUDIO | c.SDL_INIT_VIDEO | c.SDL_INIT_EVENTS) != 0) {
        c.SDL_Log("Unable to initialize SDL: %s", c.SDL_GetError());
        return error.SDLInitializationFailed;
    }

    var flags = c.SDL_WINDOW_OPENGL | c.SDL_WINDOW_RESIZABLE;
    if (config.fullscreen)
        flags |= c.SDL_WINDOW_FULLSCREEN_DESKTOP;

    window = c.SDL_CreateWindow(
        "zig-engine",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        @as(c_int, @intCast(config.windowWidth)),
        @as(c_int, @intCast(config.windowHeight)),
        @as(u32, @intCast(flags)),
    ) // c.SDL_WINDOW_FULLSCREEN_DESKTOP ) //c.SDL_WINDOW_RESIZABLE)
        orelse
        {
            c.SDL_Log("Unable to create window: %s", c.SDL_GetError());
            return error.SDLInitializationFailed;
        };

    renderer = c.SDL_CreateRenderer(
        window,
        -1,
        c.SDL_RENDERER_ACCELERATED,
    ) orelse
        {
            c.SDL_Log("Unable to create renderer: %s", c.SDL_GetError());
            return error.SDLInitializationFailed;
        };

    renderTexture = c.SDL_CreateTexture(
        renderer,
        c.SDL_PIXELFORMAT_ABGR8888,
        c.SDL_TEXTUREACCESS_STATIC,
        @as(c_int, @intCast(config.renderWidth)),
        @as(c_int, @intCast(config.renderHeight)),
    );

    uiSDLBackend = SDLBackend{
        .window = @as(*SDLBackend.c.SDL_Window, @ptrCast(window)),
        .renderer = @as(*SDLBackend.c.SDL_Renderer, @ptrCast(renderer)),
    };

    // init dvui Window (maps onto a single OS window)
    uiWindow = try dvui.Window.init(
        @src(),
        gpa,
        uiSDLBackend.backend(),
        .{},
    );
}

pub fn showMouseCursor(show: u1) u1 {
    return @as(u1, @intCast(c.SDL_ShowCursor(@as(c_int, @intCast(show)))));
}

pub fn setRelativeMouseMode(enabled: u1) void {
    _ = c.SDL_SetRelativeMouseMode(enabled);
}

pub fn setCaptureMouse(enabled: u1) void {
    _ = c.SDL_CaptureMouse(enabled);
}

pub inline fn targetFrameTimeMs() u32 {
    return config.targetDt();
}

/// Set Render texture pointer
pub fn updateRenderTexture(data: *u8, len: usize) void {
    const pixelsPtr = @as(*anyopaque, @ptrCast(data));
    if (c.SDL_UpdateTexture(renderTexture, 0, pixelsPtr, @as(c_int, @intCast(len))) != 0)
        c.SDL_Log("Unable to update texture: %s", c.SDL_GetError());
}

/// Release resources
pub fn shutdown() void {
    uiWindow.deinit();
    common.shutdown();

    c.SDL_DestroyRenderer(renderer);
    c.SDL_DestroyWindow(window);
    c.SDL_Quit();
}

/// Start updating the system
/// returns false when system quit message handled
pub fn beginUpdate() bool {
    const tracy = trace(@src());
    defer tracy.end();

    _ = common.beginUpdate();

    t0 = @as(u32, @intCast(c.SDL_GetTicks()));

    uiWindow.begin(std.time.nanoTimestamp()) catch unreachable;

    var event: c.SDL_Event = undefined;
    while (c.SDL_PollEvent(&event) != 0) {
        if (event.type == c.SDL_QUIT) {
            return false;
        }

        if (uiSDLBackend.addEvent(&uiWindow, event) catch continue)
            continue;

        switch (event.type) {
            c.SDL_QUIT => {
                return false;
            },

            c.SDL_KEYDOWN, c.SDL_KEYUP => {
                const keyCode = @as(input.KeyCode, @enumFromInt(@as(u16, @intCast(event.key.keysym.scancode))));

                input.setKeyState(keyCode, @as(u1, @intCast(event.key.state)));
            },

            c.SDL_MOUSEMOTION => {
                input.setMousePos(event.motion.x, event.motion.y);
            },

            c.SDL_MOUSEBUTTONDOWN => {
                input.setMouseButton(event.button.button, 1);
            },

            c.SDL_MOUSEBUTTONUP => {
                input.setMouseButton(event.button.button, 0);
            },

            else => {},
        }
    }

    return true;
}

/// Copy render texture to device
pub fn copyRenderTexture() void {
    _ = c.SDL_RenderClear(renderer);
    if (c.SDL_RenderCopy(renderer, renderTexture, 0, 0) != 0)
        c.SDL_Log("Unable to copy texture: %s", c.SDL_GetError());
}

/// Copy render texture to device
pub fn renderPresent() void {
    copyRenderTexture();
    _ = uiWindow.end(.{}) catch unreachable;

    _ = c.SDL_RenderPresent(renderer);
}

/// Stop processing this frame
pub fn endUpdate() u32 {
    common.endUpdate();

    const t1 = @as(u32, @intCast(c.SDL_GetTicks()));
    const dtInt = t1 - t0;
    // if (dtInt < config.targetDt())
    //     c.SDL_Delay((config.targetDt() - dtInt) - 1);

    t0 = t1;

    return dtInt;
}
