const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const run_step = b.step("run", "Run the demo");

    var exe = b.addExecutable(.{
        .name = "example",
        .root_source_file = .{ .path = "src/example.zig" },
        .target = target,
        .optimize = optimize,
    });

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());

    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    b.getInstallStep().dependOn(&b.addInstallArtifact(exe, .{ .dest_dir = .{ .override = .{ .custom = "../bin" } } }).step);

    run_step.dependOn(&run_cmd.step);

    const zpoly_mod_dep = b.dependency("zpoly_zon_name", .{ // as declared in build.zig.zon
        .target = target,
        .optimize = optimize,
    });

    const zpoly_mod = zpoly_mod_dep.module("zpoly_import_name"); // as declared in build.zig of dependency

    exe.root_module.addImport("zpoly_import_name", zpoly_mod); // name to use when importing

}
