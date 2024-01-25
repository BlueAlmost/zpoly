const std = @import("std");
const print = std.debug.print;
const zpoly = @import("zpoly_import_name");
const Complex = std.math.Complex;

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const R = f64;
    const C = Complex(R);

    var f = try zpoly.Polynomial(R).init(allocator, 8, 6);

    f.zeros();
    f.val[f.val.len - 1] = 1;
    f.val[0] = -1;

    print("\npolynomial f: \n", .{});
    print("f: {any}\n", .{f});
    f.prt();

    const x = try allocator.alloc(C, f.val.len - 1);
    try zpoly.roots(R, allocator, f, x);

    print("roots:\n", .{});
    for (x, 0..) |val, i| {
        print("\ti: {d}, ({d:.4}, {d:.4})\n", .{ i, val.re, val.im });
    }
}
