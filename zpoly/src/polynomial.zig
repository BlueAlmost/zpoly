const std = @import("std");
const print = std.debug.print;
const math = std.math;
const Complex = std.math.complex.Complex;

const zpoly = @import("../zpoly.zig");
const cmath = zpoly.cmath;

const fmtCi = zpoly.fmtCi;
const fmtCp = zpoly.fmtCp;

pub fn prt_ascend(comptime T: type, comptime fmtstr: []const u8, a: []T) void {
    // print representation as: a_n * x^n + ... + a_1 * x + a_0
    switch (T) {
        f32, f64 => {
            for (a, 0..) |val, i| {
                if (val != 0) {
                    if (i == 0) {
                        print(fmtstr, .{val});
                    } else {
                        if (val < 0) {
                            print(" - " ++ fmtstr ++ " x^{d}", .{ -val, i });
                        } else {
                            print(" + " ++ fmtstr ++ " x^{d}", .{ val, i });
                        }
                    }
                }
            }
        },

        Complex(f32), Complex(f64) => {
            for (a, 0..) |val, i| {
                if ((val.re != 0) or (val.im != 0)) {
                    if (i == 0) {
                        print(fmtstr, .{fmtCp(val)});
                    } else {
                        print(" + " ++ fmtstr ++ " x^{d}", .{ fmtCp(val), i });
                    }
                }
            }
        },

        else => @compileError("requested Polynomial type is not allowed"),
    }
    print(" \n", .{});
}

pub fn prt_descend(comptime T: type, comptime fmtstr: []const u8, a: []T) void {
    // print representation as: a_n * x^n + ... + a_1 * x + a_0
    switch (T) {
        f32, f64 => {
            var i: usize = a.len;
            if (i == 1) {
                print(fmtstr, .{a[i - 1]});
            } else {
                var first: bool = true;
                while (i != 0) {
                    i -= 1;
                    const val: T = a[i];

                    if (val != 0) {
                        if (first) {
                            print(fmtstr ++ " x^{d}", .{ val, i });
                            first = false;
                        } else {
                            if (val < 0) {
                                print(" - " ++ fmtstr, .{-val});
                                if (i != 0) print(" x^{d}", .{i});
                            } else if (val > 0) {
                                print(" + " ++ fmtstr, .{val});
                                if (i != 0) print(" x^{d}", .{i});
                            }
                        }
                    }
                }
            }
        },

        Complex(f32), Complex(f64) => {
            var i: usize = a.len;
            if (i == 1) {
                print(fmtstr, .{fmtCp(a[i - 1])});
            } else {
                var first: bool = true;
                while (i != 0) {
                    i -= 1;
                    const val: T = a[i];

                    if ((val.re != 0) or (val.im != 0)) {
                        if (first) {
                            print(fmtstr ++ " x^{d}", .{ fmtCp(val), i });
                            first = false;
                        } else {
                            print(" + " ++ fmtstr, .{fmtCp(val)});
                            if (i != 0) print(" x^{d}", .{i});
                        }
                    }
                }
            }
        },

        else => @compileError("requested Polynomial type is not allowed"),
    }
    print(" \n", .{});
}

pub fn fill(comptime T: type, x: []T, new_val: T) void {
    for (x) |*val| {
        val.* = new_val;
    }
}

pub fn monic(comptime T: type, x: []T) void {
    switch (T) {
        f32, f64 => {
            const alpha: T = 1 / x[x.len - 1];
            for (x) |*val| {
                val.* *= alpha;
            }
            x[x.len - 1] = 1; // fix value to exactly unity

        },

        Complex(f32), Complex(f64) => {
            const R = @TypeOf(x[0].re);

            const lead: T = x[x.len - 1];
            const beta: R = lead.re * lead.re + lead.im * lead.im;
            const alpha = T.init(lead.re / beta, -lead.im / beta);

            for (x) |*val| {
                val.* = cmath.mul(alpha, val.*);
            }
            // fix value to exactly unity
            x[x.len - 1] = T.init(1, 0);
        },
        else => {
            @compileError("unknown Polynomial type");
        },
    }
}

pub fn neg(comptime T: type, x: []T) void {
    switch (T) {
        f32, f64 => {
            for (x) |*val| {
                val.* = -val.*;
            }
        },

        Complex(f32), Complex(f64) => {
            for (x) |*val| {
                val.* = cmath.neg(val.*);
            }
        },
        else => @compileError("requested Polynomial type is not allowed"),
    }
}

pub fn ones(comptime T: type, x: []T) void {
    switch (T) {
        f32, f64 => {
            for (x) |*val| {
                val.* = 1;
            }
        },
        Complex(f32), Complex(f64) => {
            for (x) |*val| {
                val.* = T.init(1, 0);
            }
        },
        else => @compileError("requested Polynomial type is not allowed"),
    }
}

pub fn scale(comptime T: type, x: []T, alpha: anytype) void {
    comptime {
        if ((@TypeOf(alpha) != f32) and (@TypeOf(alpha) != f64) and (@TypeOf(alpha) != comptime_float)) {
            @compileError("scaling parameter must be real");
        }
    }

    switch (T) {
        f32, f64 => {
            for (x) |*val| {
                val.* *= alpha;
            }
        },

        Complex(f32), Complex(f64) => {
            for (x) |*val| {
                val.* = cmath.scale(val.*, alpha);
            }
        },
        else => @compileError("requested Polynomial type is not allowed"),
    }
}

pub fn zeros(comptime T: type, x: []T) void {
    switch (T) {
        f32, f64 => {
            for (x) |*val| {
                val.* = 0;
            }
        },
        Complex(f32), Complex(f64) => {
            for (x) |*val| {
                val.* = T.init(0, 0);
            }
        },
        else => @compileError("requested Polynomial type is not allowed"),
    }
}

test "\t fill - real \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(R, 3);
        fill(R, x, 12.34);

        try std.testing.expectApproxEqAbs(@as(R, 12.34), x[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 12.34), x[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, 12.34), x[2], eps);
    }
}

test "\t fill - complex \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(C, 3);
        fill(C, x, C.init(1.23, 4.56));

        try std.testing.expectApproxEqAbs(@as(R, 1.23), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.56), x[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.23), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.56), x[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.23), x[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.56), x[2].im, eps);
    }
}

test "\t monic - real \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(R, 4);
        x[0] = 2.1;
        x[1] = 4.2;
        x[2] = -0.21;
        x[3] = 2.1;

        monic(R, x);

        try std.testing.expectApproxEqAbs(@as(R, 1.0), x[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 2.0), x[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.1), x[2], eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.0), x[3], eps);
    }
}

test "\t monic - complex \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(C, 4);
        x[0].re = 2.1;
        x[0].im = 2.1;

        x[1].re = 4.2;
        x[1].im = -21.0;

        x[2].re = 0.21;
        x[2].im = 0.21;

        x[3].re = 2.1;
        x[3].im = -3.3;

        monic(C, x);

        try std.testing.expectApproxEqAbs(@as(R, -0.164705882353), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.741176405888), x[0].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 5.105882352941), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -1.976470588235), x[1].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, -0.016470588235), x[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.074117647058), x[2].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 1.0), x[3].re, eps);
    }
}

test "\t neg - real \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(R, 2);
        x[0] = 1.23;
        x[1] = 4.56;
        neg(R, x);

        try std.testing.expectApproxEqAbs(@as(R, -1.23), x[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, -4.56), x[1], eps);
    }
}

test "\t neg - complex \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(C, 2);
        x[0] = C.init(1.23, 4.56);
        x[1] = C.init(7.89, -4.56);
        neg(C, x);

        try std.testing.expectApproxEqAbs(@as(R, -1.23), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -4.56), x[0].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, -7.89), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.56), x[1].im, eps);
    }
}

test "\t ones - real \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(R, 2);
        ones(R, x);

        try std.testing.expectApproxEqAbs(@as(R, 1.0), x[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.0), x[1], eps);
    }
}

test "\t ones - complex \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(C, 2);
        ones(C, x);

        try std.testing.expectApproxEqAbs(@as(R, 1.0), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[0].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 1.0), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[1].im, eps);
    }
}

test "\t scale - real \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(R, 2);
        x[0] = 1.1;
        x[1] = 3.3;

        const alpha: R = 2.0;

        scale(R, x, alpha);

        try std.testing.expectApproxEqAbs(@as(R, 2.2), x[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 6.6), x[1], eps);
    }
}

test "\t scale - complex \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        comptime var C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(C, 2);
        x[0] = C.init(1.1, 3.3);
        x[1] = C.init(1.1, -3.3);

        // var alpha: R = 2.0;
        // x.scale(alpha);
        scale(C, x, 2.0);

        try std.testing.expectApproxEqAbs(@as(R, 2.2), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 6.6), x[0].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 2.2), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -6.6), x[1].im, eps);
    }
}

test "\t zeros - real \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(R, 3);
        zeros(R, x);

        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[2], eps);
    }
}

test "\t zeros - complex \n" {
    const eps = 1e-6;

    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const x = try allocator.alloc(C, 3);
        zeros(C, x);

        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[0].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[1].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[2].im, eps);
    }
}
