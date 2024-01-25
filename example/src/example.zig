const std = @import("std");
const print = std.debug.print;
const zpoly = @import("zpoly_import_name");
const Complex = std.math.Complex;

const fmtCp = zpoly.fmtCp;

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const R = f64;
    const C = Complex(R);

    print("\nREAL POLYNOMIALS ###########################################################\n", .{});

    // instantiating and printing real polynomials -----------------------------
    var f = try allocator.alloc(R, 4);
    f[0..4].* = .{ 1, 1, 0.0, 2 };
    print("\npolynomial f(x): \n", .{});
    zpoly.prt_descend(R, "{d:.4}", f);
    zpoly.prt_ascend(R, "{d:.4}", f);

    var g = try allocator.alloc(R, 4);
    g[0..4].* = .{ -1, -2, 3, 4 };
    print("\npolynomial g(x): \n", .{});
    zpoly.prt_descend(R, "{d:.4}", g);
    zpoly.prt_ascend(R, "{d:.4}", g);

    // // make g a monic polynomial -----------------------------------------------
    print("\npolynomial g(x) transformed to a monic polynomial: \n", .{});
    zpoly.monic(R, g);
    zpoly.prt_descend(R, "{d:.4}", g);
    zpoly.prt_ascend(R, "{d:.4}", g);

    // // add two polynomials -----------------------------------------------------
    const a = try allocator.alloc(R, 4);
    try zpoly.add(R, a, f, g);
    print("\npolynomial a(x) = f(x) + g(x): \n", .{});
    zpoly.prt_descend(R, "{d:.4}", a);
    zpoly.prt_ascend(R, "{d:.4}", a);

    // // multiply two polynomials ------------------------------------------------
    const b = try allocator.alloc(R, 7);
    try zpoly.mult(R, b, f, g);
    print("\npolynomial b(x) = f(x) * g(x): \n", .{});
    zpoly.prt_descend(R, "{d:.4}", b);
    zpoly.prt_ascend(R, "{d:.4}", b);

    // // evaluate a polynomial at a point ----------------------------------------
    print("\nevaluate polynomial b(x) at x=2: {d}\n", .{zpoly.eval(R, R, b, 2)});

    // determine all roots of a polynomial --------------------------------------
    const x = try allocator.alloc(C, b.len - 1); // alloc memory for root values
    zpoly.prt_descend(R, "{d:.4}", b);
    try zpoly.roots(R, allocator, b, x);
    print("\nroots of polynomial b(x):\n", .{});
    for (x, 0..) |val, i| {
        print("\ti: {d}, {d:.4})\n", .{ i, fmtCp(val) });
    }

    print("\nevaluate polynomial at its roots:\n", .{});
    for (x, 0..) |val, i| {
        const b_eval = zpoly.eval(R, C, b, val);
        print("\tvalue at root i: {d}, {d:.4}\n", .{ i, fmtCp(b_eval) });
    }

    // set lower order coefficients to zero (to verify proper root finding when
    // there are multiple roots at the origin
    b[0] = 0;
    b[1] = 0;
    print("\n", .{});
    zpoly.prt_descend(R, "{d:.4}", b);
    try zpoly.roots(R, allocator, b, x);
    print("\nroots of polynomial b(x):\n", .{});
    for (x, 0..) |val, i| {
        print("\ti: {d}, {d:.4})\n", .{ i, fmtCp(val) });
    }
    print("\nevaluate polynomial at its roots:\n", .{});
    for (x, 0..) |val, i| {
        const b_eval = zpoly.eval(R, C, b, val);
        print("\tvalue at root i: {d}, {d:.4}\n", .{ i, fmtCp(b_eval) });
    }

    print("\nCOMPLEX POLYNOMIALS #######################################################\n", .{});

    // instantiating and printing complex polynomials --------------------------
    var cf = try allocator.alloc(C, 4);
    cf[0..4].* = .{ C.init(1, 2), C.init(-1, 3), C.init(0, 0), C.init(2, -1) };
    print("\npolynomial cf(x): \n", .{});
    zpoly.prt_descend(C, "{d:.4}", cf);
    zpoly.prt_ascend(C, "{d:.4}", cf);

    var cg = try allocator.alloc(C, 4);
    cg[0..4].* = .{ C.init(-1, 3), C.init(-2, 1), C.init(3, 2), C.init(4, 1) };
    print("\npolynomial cg(x): \n", .{});
    zpoly.prt_descend(C, "{d:.4}", cg);
    zpoly.prt_ascend(C, "{d:.4}", cg);

    // // make cg a monic polynomial -----------------------------------------------
    print("\npolynomial cg(x) transformed to a monic polynomial: \n", .{});
    zpoly.monic(C, cg);
    zpoly.prt_descend(C, "{d:.4}", cg);
    zpoly.prt_ascend(C, "{d:.4}", cg);

    // // add two polynomials -----------------------------------------------------
    const ca = try allocator.alloc(C, 4);
    try zpoly.add(C, ca, cf, cg);
    print("\npolynomial ca(x) = cf(x) + cg(x): \n", .{});
    zpoly.prt_descend(C, "{d:.4}", cg);
    zpoly.prt_ascend(C, "{d:.4}", cg);

    // // multiply two polynomials ------------------------------------------------
    const cb = try allocator.alloc(C, 7);
    try zpoly.mult(C, cb, cf, cg);
    print("\npolynomial cb(x) = cf(x) * cg(x): \n", .{});
    zpoly.prt_descend(C, "{d:.4}", cb);
    zpoly.prt_ascend(C, "{d:.4}", cb);

    // evaluate a complex polynomial at a real point ----------------------------
    const rx = 2;
    const cb2r = zpoly.eval(C, R, cb, rx);
    print("\nevaluate complex polynomial at real value, cb({d:.4}): {d:.4}\n", .{ rx, fmtCp(cb2r) });

    // evaluate a complex polynomial at a complex point ----------------------------
    const cx = C.init(1.1, -2.2);
    const cb2c = zpoly.eval(C, C, cb, cx);
    print("\nevaluate complex polynomial at complex value, cb({d:.4}): {d:.4}\n", .{ fmtCp(cx), fmtCp(cb2c) });

    // determine all roots of a complex polynomial --------------------------------------
    const cx_roots = try allocator.alloc(C, cb.len - 1); // alloc memory for root values
    zpoly.prt_descend(C, "{d:.4}", cb);
    try zpoly.roots(C, allocator, cb, cx_roots);
    print("\nroots of polynomial cb:\n", .{});
    for (cx_roots, 0..) |val, i| {
        print("\ti: {d}, {d:.4}\n", .{ i, fmtCp(val) });
    }

    print("\nevaluate complex polynomical at its roots:\n", .{});
    for (cx_roots, 0..) |val, i| {
        const cb_eval = zpoly.eval(C, C, cb, val);
        print("\tvalue at root i: {d}, {d:.4}\n", .{ i, fmtCp(cb_eval) });
    }

    // set lower order coefficients to zero (to verify proper root finding when
    // there are multiple roots at the origin
    cb[0] = C.init(0, 0);
    cb[1] = C.init(0, 0);
    print("\n", .{});
    zpoly.prt_descend(C, "{d:.4}", cb);

    try zpoly.roots(C, allocator, cb, cx_roots);
    print("\nroots of polynomial cb:\n", .{});
    for (cx_roots, 0..) |val, i| {
        print("\ti: {d}, {d:.4}\n", .{ i, fmtCp(val) });
    }

    print("\nevaluate complex polynomical at its roots:\n", .{});
    for (cx_roots, 0..) |val, i| {
        const cb_eval = zpoly.eval(C, C, cb, val);
        print("\tvalue at root i: {d}, {d:.4}\n", .{ i, fmtCp(cb_eval) });
    }
}
