const std = @import("std");
const warn = std.debug.print;
const fmt = std.fmt;
const assert = @import("std").debug.assert;
const math = std.math;
const io = std.io;
const Allocator = std.mem.Allocator;

const Vec4f = @import("../../core/vector.zig").Vec4f;
const Mat44f = @import("../../core/matrix.zig").Mat44f;
const Profile = @import("../../core/profiler.zig").Profile;

const Mesh = @import("../../render/mesh.zig").Mesh;

const NewLine = 0x0A;
const WhiteSpace = 0x20;
const ForwardSlash = '/';

pub const State = enum {
    // These must be first with these explicit values as we rely on them for indexing the
    // bit-stack directly and avoiding a branch.
    Read,
};

pub fn find(buf: []const u8, start: usize, char: u8) ![]const u8 {
    var i: usize = start;
    var found = false;
    while (i < buf.len) {
        const c = buf[i];
        if (c == char) {
            found = true;
        } else if (found) {
            return buf[start .. i - 1];
        }

        i += 1;
    }

    return buf;
}

pub fn parseVertex(buf: []const u8) ![2]Vec4f {
    var iter = std.mem.split(u8, buf, " ");
    var count: usize = 0;
    var pos = Vec4f.init(0, 0, 0, 1);
    var color = Vec4f.one();
    while (iter.next()) |el| {
        if (el.len < 1)
            continue;

        count += 1;

        switch (count) {
            // position
            1 => {
                pos.setX(try std.fmt.parseFloat(f32, el));
            },
            2 => {
                pos.setY(try std.fmt.parseFloat(f32, el));
            },
            3 => {
                pos.setZ(try std.fmt.parseFloat(f32, el));
            },

            // color
            4 => {
                color.setX(try std.fmt.parseFloat(f32, el));
            },
            5 => {
                color.setY(try std.fmt.parseFloat(f32, el));
            },
            6 => {
                color.setZ(try std.fmt.parseFloat(f32, el));
            },
            7 => {
                color.setW(try std.fmt.parseFloat(f32, el));
            },
            else => {},
        }
    }

    return .{ pos, color };
}

pub fn parseVector(buf: []const u8) !Vec4f {
    var iter = std.mem.split(u8, buf, " ");
    var count: usize = 0;
    var vec = Vec4f.zero();
    while (iter.next()) |el| {
        if (el.len < 1)
            continue;

        count += 1;

        switch (count) {
            // position
            1 => {
                vec.setX(try std.fmt.parseFloat(f32, el));
            },
            2 => {
                vec.setY(try std.fmt.parseFloat(f32, el));
            },
            3 => {
                vec.setZ(try std.fmt.parseFloat(f32, el));
            },
            4 => {
                vec.setW(try std.fmt.parseFloat(f32, el));
            },
            else => {},
        }
    }

    return vec;
}

const FacePointDef = struct {
    v: u16,
    uv: u16,
    vn: u16,
};

pub fn parseFacePoint(buf: []const u8) !FacePointDef {
    // var data = [_]u16{0} ** 3;

    // vertex index
    const vstr = try find(buf, 0, ForwardSlash);
    var nstart = vstr.len + 1;

    // uv coord index
    const uvstr = try find(buf, vstr.len + 1, ForwardSlash);

    // vertex normal index
    nstart += uvstr.len + 1;
    const nstr = buf[nstart..buf.len];

    return FacePointDef{
        .v = try std.fmt.parseInt(u16, vstr, 10),
        .uv = try std.fmt.parseInt(u16, uvstr, 10),
        .vn = try std.fmt.parseInt(u16, nstr, 10),
    };
}

pub fn parseFace(buf: []const u8) !struct {
    v0: FacePointDef,
    v1: FacePointDef,
    v2: FacePointDef,
} {
    // var data = [_]u16{0} ** 9;

    const xstr = try find(buf, 0, WhiteSpace);
    const ystr = try find(buf, xstr.len + 1, WhiteSpace);
    const zstart = xstr.len + 1 + ystr.len + 1;
    const zstr = buf[zstart..buf.len];

    return .{
        .v0 = (try parseFacePoint(xstr)),
        .v1 = (try parseFacePoint(ystr)),
        .v2 = (try parseFacePoint(zstr)),
    };
}

pub fn importObjFile(allocator: *Allocator, file_path: []const u8) !Mesh {
    const cwd = std.fs.cwd();

    const resolvedPath = try std.fs.path.resolve(allocator.*, &[_][]const u8{file_path});
    defer allocator.free(resolvedPath);

    std.debug.print("obj path: {s}\n", .{resolvedPath});

    var file = try cwd.openFile(resolvedPath, .{});
    defer file.close();

    var stream_source = io.StreamSource{ .file = file };
    var in = stream_source.reader();
    // var state: State = .Read;
    // var data: u8 = 0;

    var buf = [_]u8{0} ** 128;

    var verts = std.ArrayList(Vec4f).init(allocator.*);
    defer verts.deinit();

    var colors = std.ArrayList(Vec4f).init(allocator.*);
    defer colors.deinit();

    var normals = std.ArrayList(Vec4f).init(allocator.*);
    defer normals.deinit();

    var texCoords = std.ArrayList(f32).init(allocator.*);
    defer texCoords.deinit();

    var triVerts = std.ArrayList(u16).init(allocator.*);
    defer triVerts.deinit();

    var triUvs = std.ArrayList(u16).init(allocator.*);
    defer triUvs.deinit();

    var triNormals = std.ArrayList(u16).init(allocator.*);
    defer triNormals.deinit();
    var lineCount: usize = 0;
    var eof = false;
    while (true) {
        var read: usize = 0;
        var line: []const u8 = undefined;
        while (true) {
            // read next line
            const byte = in.readByte() catch |err| switch (err) {
                error.EndOfStream => {
                    eof = true;
                    break;
                },
                else => {
                    eof = true;
                    break;
                },
            };

            if (byte == NewLine) {
                // also skip return character if it's there
                if (read > 0 and buf[read - 1] == '\r')
                    read -= 1;

                line = buf[0..read];
                lineCount += 1;
                break;
            }

            buf[read] = byte;
            read += 1;
        }

        if (eof and read == 0)
            break;

        if (read <= 4)
            continue;

        std.debug.print("{}: {s}\n", .{ lineCount, line });
        // std.debug.warn("[{any}]: {any}\n", .{lineCount, line});
        switch (line[0]) {
            // vertex information
            'v' => {
                switch (line[1]) {
                    // vertex
                    WhiteSpace => {
                        const parsed = try parseVertex(line[2..]);
                        try verts.append(parsed[0]);
                        try colors.append(parsed[1]);
                    },

                    // texture coordinate
                    't' => {
                        const uvs = try parseVector(line[3..]);
                        try texCoords.append(uvs.x());
                        try texCoords.append(uvs.y());
                    },

                    // vertex normal
                    'n' => {
                        const normal = (try parseVector(line[3..])).normalized3();
                        try normals.append(normal);
                    },
                    else => continue,
                }
            },

            // face information
            'f' => {
                const faceData = try parseFace(line[2..]);

                try triVerts.append(faceData.v0.v - 1);
                try triUvs.append(faceData.v0.uv - 1);
                try triNormals.append(faceData.v0.vn - 1);

                try triVerts.append(faceData.v1.v - 1);
                try triUvs.append(faceData.v1.uv - 1);
                try triNormals.append(faceData.v1.vn - 1);

                try triVerts.append(faceData.v2.v - 1);
                try triUvs.append(faceData.v2.uv - 1);
                try triNormals.append(faceData.v2.vn - 1);
            },

            // Object name
            'o' => {},

            's' => {},

            // usemtl <material name>
            'u' => {},

            // mtllib <material file path>
            'm' => {},

            else => {
                // skip this line
                continue;
            },
        }
    }

    return Mesh.init(
        try verts.toOwnedSlice(),
        try triVerts.toOwnedSlice(),
        try triNormals.toOwnedSlice(),
        try triUvs.toOwnedSlice(),
        try colors.toOwnedSlice(), // vertex colors
        try normals.toOwnedSlice(),
        try texCoords.toOwnedSlice(),
    );
}

// test "Mesh .obj file load"
// {
//   // importObjFile()
// }
