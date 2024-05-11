const std = @import("std");
const vector = @import("vector.zig");
const Vec4f = vector.Vec4f;
const assert = std.debug.assert;
const StructField = std.builtin.TypeInfo.StructField;

pub const Register = struct { 

  pub const Definition = struct {
    vt: type, // value type
    vd: anytype, // value default
    rt: Runtime.Flags,  // should sample the value
    info: Info,
  };

  pub const Info = struct {
    pub const Flags = enum(u8) {
      Read = 1 << 0,
      Write = 1 << 1,
      Save = 1 << 2,
      Load = 1 << 3,
      
      
      // const RW = Read | Write;
      // const Config = Save | Load;

    };
    
    description:[] const u8,
    flags:Flags=Flags.Read,
    meta: anytype = {}
  };

  pub const Runtime = struct {
    pub const Flags = enum(u8) {
      None = 0,
      Sample = 1 << 0,
    };
  };

  pub fn define(comptime valueType:type, valueDefault:valueType, runtimeFlags:Runtime.Flags, info:Info) Definition {
    return Definition{ .vt = valueType, .vd = valueDefault, .rt = runtimeFlags, .info = info };
  }

  pub fn Type(
    name: []const u8, 
    comptime T: type, 
    default: T, 
    runtimeFlags:Runtime.Flags, 
    info: Info
  ) type {
    return struct {
        const Self = @This();
        pub const info = info;

        value: T = default,
        flags: Runtime.Flags = runtimeFlags,

        pub inline fn name(self: Self) []const u8 {
            _ = self;
            return name;
        }

        pub inline fn default(self: Self) T {
            _ = self;
            return default;
        }

        pub fn get(self: Self) T {
            return self.value;
        }

        pub fn set(self: *Self, v: T) void {
            self.value = v;
        }

        pub fn reset(self: *Self) void {
            self.value = self.default;
        }
    };
  }  
};


const categoryName = "isCategory";

fn buildRegisters(comptime T: anytype, comptime isCategory: bool) []const StructField {    
    comptime var fields: []const StructField = &[_]StructField{};

    inline for (@typeInfo(T).Struct.fields) |field| {
        if (field.field_type == Register.Definition) 
        {
            const value = field.default_value.?;
            const regType = Register.Type(field.name, value.vt, value.vd, value.rt, value.info);
            const regTypeDefault = regType{};

            fields = fields ++ &[_]StructField{.{ 
              .name = field.name, 
              .field_type = regType, 
              .default_value = regTypeDefault, 
              .is_comptime = false, 
              .alignment = @alignOf(field.field_type) 
              }};
        } 
        else 
        {
            // * Build a new Category
            const Category = @Type(.{ .Struct = .{
                .layout = .Auto,
                .fields = buildRegisters(field.default_value.?, true),
                .decls = &[_]std.builtin.TypeInfo.Declaration{},
                .is_tuple = false,
            } });

            const CategoryDefault = Category{};

            fields = fields ++ &[_]StructField{.{ 
              .name = field.name, 
              .field_type = Category, 
              .default_value = CategoryDefault, 
              .is_comptime = false, 
              .alignment = @alignOf(T) 
            }};
        }
    }

    // Add isCateogry flag
    fields = fields ++ &[_]StructField{.{ 
      .name =categoryName, 
      .field_type = bool, 
      .default_value = isCategory, 
      .is_comptime = false, 
      .alignment = @alignOf(T) 
      }
    };

    return fields;
}



pub fn Telemetry(comptime T: anytype) type {

    // create a new type based on all of the fields of T
    // each T field should be a ValueDefinition

    comptime var fields = buildRegisters(T, true);

    comptime var RegistersType = @Type(.{ .Struct = .{
        .layout = .Auto,
        .fields = fields,
        .decls = &[_]std.builtin.TypeInfo.Declaration{},
        .is_tuple = false,
    } });

    // @compileLog("fields", @typeInfo(T).Struct.decls.len, fields.len);

    return struct {
        const Self = @This();

        // all Value types
        reg: RegistersType = .{},

        pub fn init() Self {
            return Self{};
        }

        pub fn printValue(
          comptime parentName: []const u8, 
          comptime name: []const u8, 
          data: anytype, 
          comptime level: u8, 
          comptime levelStep:u8) !void {
            const stdout = std.io.getStdOut().writer();
            const hasParentName = parentName.len > 0;
            const indentChar = " ";

            if (@hasField(@TypeOf(data), categoryName)) {
                if (hasParentName)
                    try stdout.print("{s}[{s}]\n", .{ indentChar ** level, name });

                inline for (std.meta.fields(@TypeOf(data))) |field| {
                    if (comptime !std.mem.eql(u8, field.name, categoryName)) {
                        const prefix = if (hasParentName) parentName ++ "." ++ field.name else field.name;
                        try printValue(prefix, field.name, @field(data, field.name), level + levelStep, levelStep);
                    }
                }
            } else {
                try stdout.print("{s}{s} = {}\n", .{ indentChar ** level, parentName, data.value });
            }
        }

        pub fn print(self: Self, comptime indent:u8) !void {
            _ = self;


            try printValue("", "", self.reg, 0, indent);
        }
    };
}

test "Telemetry 1" {
    std.debug.print("\n", .{});

    const rtf = Register.Runtime.Flags;
    const info = Register.Info.Flags;
    const types2 = struct {
        motor: anytype = struct { 
          output: anytype = struct {
            enabled:Register.Definition = Register.define(bool, false, rtf.None, .{ 
              .flags= info.Write, .description = "sadlfsf" }),} 
        },
        fc: anytype = struct {
            armed: Register.Definition = Register.define(bool, false, rtf.None, .{ .description = "sadlfsf" } ),
        },
        dev: anytype = struct { 
          imu: anytype = struct { 
            gyro: anytype = struct {
              rates: Register.Definition = Register.define(Vec4f, Vec4f.zero, rtf.Sample, .{ .description = "Gyro rates deg/s" } ),
        } } },
    };

    var t = Telemetry(types2).init();

    try t.print(1);

    assert(t.reg.dev.imu.gyro.rates.value.x == 0.0);

    var gyro = t.reg.dev.imu.gyro;
    gyro.rates.value.x = 1.0;
    // t.reg.dev.imu.gyro.rates.value.x = 1.0;

    assert(t.reg.dev.imu.gyro.rates.value.x == 1.0);
    try t.print(4);
}
