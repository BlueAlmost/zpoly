const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const zpoly = b.addModule("zpoly_import_name", .{
        .root_source_file = b.path("./zpoly.zig"),
    });

    const zpoly_tests = b.addTest(.{
        .root_source_file = b.path ("./zpoly_test.zig"),
        .target = target,
        .optimize = optimize,
    });

    zpoly_tests.root_module.addImport("zpoly_import_name", zpoly);

    const run_zpoly_tests = b.addRunArtifact(zpoly_tests);

    const test_step = b.step("test", "Run module tests");

    test_step.dependOn(&run_zpoly_tests.step);
}
