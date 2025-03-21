const std = @import("std");
const builtin = @import("builtin");

// Although this function looks imperative, note that its job is to
// declaratively construct a build graph that will be executed by an external
// runner.
pub fn build(b: *std.Build) void {
    const tracy = b.option(bool, "tracy", "Enable Tracy integration. Supply path to Tracy source") orelse false;
    const tracy_callstack = b.option(bool, "tracy-callstack", "Include callstack information with Tracy data. Does nothing if -Dtracy is not provided") orelse false;
    const tracy_allocation = b.option(bool, "tracy-allocation", "Include allocation information with Tracy data. Does nothing if -Dtracy is not provided") orelse false;

    // Standard target options allows the person running `zig build` to choose
    // what target to build for. Here we do not override the defaults, which
    // means any target is allowed, and the default is native. Other options
    // for restricting supported target set are available.
    const target = b.standardTargetOptions(.{});

    // Standard optimization options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall. Here we do not
    // set a preferred release mode, allowing the user to decide how to optimize.
    const optimize = b.standardOptimizeOption(.{});

    const exe = b.addExecutable(.{
        .name = "zig-engine",
        // In this case the main source file is merely a path, however, in more
        // complicated build scripts, this could be a generated file.
        .root_source_file = b.path("src/engine.zig"),
        .target = target,
        .optimize = optimize,
    });

    const exe_options = b.addOptions();
    exe.root_module.addOptions("build_options", exe_options);
    exe_options.addOption(bool, "enable_tracy", tracy);
    exe_options.addOption(bool, "enable_tracy_callstack", tracy_callstack);
    exe_options.addOption(bool, "enable_tracy_allocation", tracy_allocation);
    exe.root_module.strip = false;

    if (tracy) {
        const tracyPath = "external/tracy/public";

        const client_cpp = std.fs.path.join(b.allocator, &[_][]const u8{ tracyPath, "TracyClient.cpp" }) catch unreachable;
        // On mingw, we need to opt into windows 7+ to get some features required by tracy.
        const tracy_c_flags: []const []const u8 = if (target.result.os.tag == .windows and target.result.abi.isGnu())
            &[_][]const u8{ "-DTRACY_ENABLE=1", "-fno-sanitize=undefined", "-D_WIN32_WINNT=0x601" }
        else
            &[_][]const u8{ "-DTRACY_ENABLE=1", "-fno-sanitize=undefined" };

        exe.addIncludePath(b.path("external/tracy/public"));
        exe.addCSourceFile(.{ .file = b.path(client_cpp), .flags = tracy_c_flags });
        exe.linkLibCpp();
        exe.linkLibC();

        if (target.result.os.tag == .windows) {
            exe.linkSystemLibrary("dbghelp");
            exe.linkSystemLibrary("ws2_32");
        }
    }

    if (builtin.os.tag == .windows) {

        //std.debug.print("b.build_root={s}",.{b.build_root.path});
        const sdl_path = b.fmt("external/win/SDL2", .{});
        // std.debug.print("{s}", .{sdl_path});
        exe.addLibraryPath(b.path(b.fmt("{s}/lib/x64", .{sdl_path})));
        exe.addIncludePath(b.path(b.fmt("{s}/include", .{sdl_path})));
        exe.addIncludePath(.{ .cwd_relative = "C:\\Program Files (x86)\\Windows Kits\\10\\Include\\10.0.17763.0\\shared\\evntprov.h" });
        b.installBinFile(b.fmt("{s}/lib/x64/SDL2.dll", .{sdl_path}), "SDL2.dll");
        exe.linkSystemLibrary("sdl2");
    } else {
        exe.linkSystemLibrary("SDL2");
    }

    const dvui_dep = b.dependency("dvui", .{ .target = target, .optimize = optimize, .backend = .sdl });

    exe.root_module.addImport("dvui", dvui_dep.module("dvui_sdl"));
    // exe.root_module.addImport("dvui", dvui_dep.module("dvui"));
    // exe.root_module.addImport("SDLBackend", dvui_dep.module("SDLBackend"));

    //exe.linkLibC();

    exe.linkSystemLibrary("c");

    // This declares intent for the executable to be installed into the
    // standard location when the user invokes the "install" step (the default
    // step when running `zig build`).
    b.installArtifact(exe);

    // This *creates* a RunStep in the build graph, to be executed when another
    // step is evaluated that depends on it. The next line below will establish
    // such a dependency.
    const run_cmd = b.addRunArtifact(exe);

    // By making the run step depend on the install step, it will be run from the
    // installation directory rather than directly from within the cache directory.
    // This is not necessary, however, if the application depends on other installed
    // files, this ensures they will be present and in the expected location.
    run_cmd.step.dependOn(b.getInstallStep());

    // This allows the user to pass arguments to the application in the build
    // command itself, like this: `zig build run -- arg1 arg2 etc`
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    // This creates a build step. It will be visible in the `zig build --help` menu,
    // and can be selected like this: `zig build run`
    // This will evaluate the `run` step rather than the default, which is "install".
    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    // Creates a step for unit testing. This only builds the test executable
    // but does not run it.
    const unit_tests = b.addTest(.{
        .root_source_file = b.path("src/engine.zig"),
        .target = target,
        .optimize = optimize,
    });

    const run_unit_tests = b.addRunArtifact(unit_tests);

    // Similar to creating the run step earlier, this exposes a `test` step to
    // the `zig build --help` menu, providing a way for the user to request
    // running the unit tests.
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_unit_tests.step);
}
