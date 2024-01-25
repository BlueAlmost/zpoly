const std = @import("std");
const print = std.debug.print;
const math = std.math;
const Complex = std.math.complex.Complex;

const zpoly = @import("../zpoly.zig");
const Polynomial = zpoly.Polynomial;
const cmath = zpoly.cmath;

const EitherCmpx = zpoly.EitherCmpx;
const ValueType = zpoly.ValueType;

pub fn fabs(comptime T: type, x: T) ValueType(T) {
    switch (T) {
        f32, f64 => {
            return @abs(x);
        },
        Complex(f32), Complex(f64) => {
            const R: type = ValueType(T);
            return std.math.hypot(R, x.re, x.im);
        },
        else => {
            @compileError("type not implemented");
        },
    }
}

pub fn add(comptime T: type, a: []T, f: []T, g: []T) !void {
    //
    // add arrays f and g of potentially different lengths.
    //
    // (useful when arrays are used to represent polynomials, i.e. have different orders).
    //    a(x) = f(x) + g(x)
    //
    //    a.len == max( f.len, g.len).

    std.debug.assert(a.len == @max(f.len, g.len));

    const min_len: usize = @min(f.len, g.len);
    const max_len: usize = @max(f.len, g.len);

    if (a.len != max_len) {
        return error.Insufficient_Len_to_Store_Result;
    }

    switch (T) {
        f32, f64 => {
            var i: usize = 0;
            while (i < min_len) : (i += 1) {
                a[i] = f[i] + g[i];
            }

            if (f.len > g.len) {
                while (i < f.len) : (i += 1) {
                    a[i] = f[i];
                }
            } else if (g.len > f.len) {
                while (i < g.len) : (i += 1) {
                    a[i] = g[i];
                }
            }
        },
        Complex(f32), Complex(f64) => {
            var i: usize = 0;
            while (i < min_len) : (i += 1) {
                a[i] = cmath.add(f[i], g[i]);
            }

            if (f.len > g.len) {
                while (i < f.len) : (i += 1) {
                    a[i] = f[i];
                }
            } else if (g.len > f.len) {
                while (i < g.len) : (i += 1) {
                    a[i] = g[i];
                }
            }
        },
        else => {
            @compileError("unsupported type");
        },
    }
}

pub fn conv(comptime T: type, a: []T, f: []T, g: []T) !void {

    // convolution of two arrays.
    //
    // useful for multiplying two polynomials (same as convolution of
    // the coefficients:  a(x) = f(x) g(x)
    //
    // a.len = f.len + g.len -1

    std.debug.assert(a.len == f.len + g.len - 1);

    if (a.len != f.len + g.len - 1) {
        return error.Insufficient_Space_to_Store_Result;
    }

    var n: usize = 0;

    switch (T) {
        f32, f64 => {
            while (n < a.len) : (n += 1) {
                a[n] = 0;
                var k: usize = 0;
                while (k < f.len) : (k += 1) {
                    if ((n >= k) and (n - k < g.len)) {
                        a[n] += f[k] * g[n - k];
                    }
                }
            }
        },
        Complex(f32), Complex(f64) => {
            while (n < a.len) : (n += 1) {
                a[n] = T.init(0, 0);

                var k: usize = 0;

                while (k < f.len) : (k += 1) {
                    if ((n >= k) and (n - k < g.len)) {
                        const tmp: T = cmath.mul(f[k], g[n - k]);
                        a[n] = cmath.add(a[n], tmp);
                    }
                }
            }
        },

        else => {
            @compileError("unsupported type");
        },
    }
}

pub fn evalDeriv(comptime A: type, comptime X: type, a: []A, x: X) EitherCmpx(A, X) {
    const R = ValueType(A);
    var n_tmp: R = undefined;

    const G = EitherCmpx(A, X);
    var g: G = undefined;

    // Evaluate using Horner's Method, the derivative of the poly whose
    // original coef array is a = [a_0 a_1 ... a_n]
    //
    // for the given value "x" as
    //
    //      g(x) = a_1 + 2*a_2*x + ... n* a_n*x^(n-1)

    switch (A) {
        f32, f64 => {
            switch (X) {
                f32, f64 => {
                    g = 0;
                    var i: usize = 1;
                    while (i < a.len) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.len - i)));
                        g = n_tmp * a[a.len - i] + x * g;
                    }
                    return g;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 1;
                    var g_tmp: G = undefined;
                    g = G.init(0, 0);
                    while (i < a.len) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.len - i)));
                        const a_tmp: A = a[a.len - i];
                        g_tmp = g;
                        g.re = n_tmp * a_tmp + (x.re * g_tmp.re - x.im * g_tmp.im);
                        g.im = (x.im * g_tmp.re + x.re * g_tmp.im);
                    }
                    return g;
                },
                else => {
                    @compileError("inconsistent polynomial and variable type");
                },
            }
        },
        Complex(f32), Complex(f64) => {
            switch (X) {
                f32, f64 => {
                    var i: usize = 1;
                    var g_tmp: G = undefined;
                    var a_tmp: G = undefined;
                    g = G.init(0, 0);
                    while (i < a.len) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.len - i)));

                        a_tmp = a[a.len - i];
                        g_tmp = g;

                        g.re = n_tmp * a_tmp.re + (x * g_tmp.re);
                        g.im = n_tmp * a_tmp.im + (x * g_tmp.im);
                    }
                    return g;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 1;
                    g = G.init(0, 0);

                    while (i < a.len) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.len - i)));

                        const tmp = cmath.mul(x, g);
                        g.re = n_tmp * a[a.len - i].re + tmp.re;
                        g.im = n_tmp * a[a.len - i].im + tmp.im;
                    }
                    return g;
                },
                else => {
                    @compileError("inconsistent polynomial and variable type");
                },
            }
        },
        else => {
            @compileError("inconsistent polynomial and variable type");
        },
    }
}

pub fn evalDerivRev(comptime A: type, comptime X: type, a: []A, x: X) EitherCmpx(A, X) {
    const R = ValueType(A);
    var n_tmp: R = undefined;

    const G = EitherCmpx(A, X);
    var g: G = undefined;

    // Evaluate using Horner's Method, the derivative of the REVERSE of the
    // poly whose original coef array is a = [a_0 a_1 ... a_n],
    //
    // thus the reversal of a is [a_n a_{n-1} ... a_0 ],
    //
    // evaluates its derivaties for the given value "x" as
    //
    //      g(x) = a_{n-1} + 2*a_{n-2}*x + ... n* a_0*x^(n-1)

    switch (A) {
        f32, f64 => {
            switch (X) {
                f32, f64 => {
                    g = 0;
                    var i: usize = 0;
                    while (i < a.len - 1) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.len - i - 1)));
                        g = n_tmp * a[i] + x * g;
                    }
                    return g;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    var g_tmp: G = undefined;
                    g = G.init(0, 0);

                    while (i < a.len - 1) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.len - i - 1)));

                        const a_tmp: A = a[i];
                        g_tmp.re = g.re;
                        g_tmp.im = g.im;

                        g.re = n_tmp * a_tmp + (x.re * g_tmp.re - x.im * g_tmp.im);
                        g.im = (x.im * g_tmp.re + x.re * g_tmp.im);
                    }
                    return g;
                },
                else => {
                    @compileError("inconsistent polynomial and variable type");
                },
            }
        },
        Complex(f32), Complex(f64) => {
            switch (X) {
                f32, f64 => {
                    var i: usize = 0;
                    g = G.init(0, 0);
                    while (i < a.len - 1) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.len - i - 1)));

                        g.re = n_tmp * a[i].re + (x * g.re);
                        g.im = n_tmp * a[i].im + (x * g.im);
                    }
                    return g;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    var g_tmp: G = undefined;
                    g = G.init(0, 0);
                    while (i < a.len - 1) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.len - i - 1)));

                        g_tmp.re = g.re;
                        g_tmp.im = g.im;

                        g.re = n_tmp * a[i].re + (x.re * g_tmp.re - x.im * g_tmp.im);
                        g.im = n_tmp * a[i].im + (x.im * g_tmp.re + x.re * g_tmp.im);
                    }
                    return g;
                },
                else => {
                    @compileError("inconsistent polynomial and variable type");
                },
            }
        },
        else => {
            @compileError("inconsistent polynomial and variable type");
        },
    }
}

pub fn eval(comptime A: type, comptime X: type, a: []A, x: X) EitherCmpx(A, X) {
    const F = EitherCmpx(A, X);
    var f: F = undefined;

    // Evaluate using Horner's Method, the polynomial having coefficent
    // array a = [a_0 a_1 ... a_n] for the given value "x" as
    //
    //      f(x) = a_0 + a_1*x + a_2*x^2 + ... a_n*x^n
    //
    // ("x" and coefficients of a must share same type)

    switch (A) {
        f32, f64 => {
            switch (X) {
                f32, f64 => {
                    f = 0;
                    var i: usize = 0;
                    while (i < a.len) : (i += 1) {
                        f = a[a.len - 1 - i] + x * f;
                    }
                    return f;
                },
                Complex(f32), Complex(f64) => {
                    var a_tmp: A = undefined;
                    var f_tmp: F = undefined;
                    f.re = 0;
                    f.im = 0;
                    var i: usize = 0;
                    while (i < a.len) : (i += 1) {
                        a_tmp = a[a.len - 1 - i];
                        f_tmp.re = f.re;
                        f_tmp.im = f.im;
                        f.re = a_tmp + (x.re * f_tmp.re - x.im * f_tmp.im);
                        f.im = (x.im * f_tmp.re + x.re * f_tmp.im);
                    }
                    return f;
                },
                else => {
                    @compileError("inconsistent polynomial and variable type");
                },
            }
        },

        Complex(f32), Complex(f64) => {
            f = F.init(0, 0);
            switch (X) {
                f32, f64 => {
                    var i: usize = 0;
                    while (i < a.len) : (i += 1) {
                        f = cmath.add(a[a.len - 1 - i], cmath.scale(f, x));
                    }
                    return f;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    while (i < a.len) : (i += 1) {
                        f = cmath.add(a[a.len - 1 - i], cmath.mul(f, x));
                    }
                    return f;
                },
                else => {
                    @compileError("inconsistent polynomial and variable type");
                },
            }
        },
        else => {
            @compileError("inconsistent polynomial and variable type");
        },
    }
}

pub fn evalRev(comptime A: type, comptime X: type, a: []A, x: X) EitherCmpx(A, X) {
    const F = EitherCmpx(A, X);
    var f: F = undefined;

    // Evaluate using Horner's Method, the REVERSED (flipped) polynomial
    // having coefficent array a = [a_0 a_1 ... a_n] for the given value "x" as
    //
    //      f(x) = a_n + a_{n-1}*x + a_{n-2}*x^2 + ... a_0*x^n
    //
    // ("x" and coefficients of a must share same type)

    switch (A) {
        f32, f64 => {
            switch (X) {
                f32, f64 => {
                    f = 0;
                    var i: usize = 0;
                    while (i < a.len) : (i += 1) {
                        f = a[i] + x * f;
                    }
                    return f;
                },
                Complex(f32), Complex(f64) => {
                    f = X.init(0, 0);
                    var i: usize = 0;
                    while (i < a.len) : (i += 1) {
                        f = cmath.mul(f, x);
                        f.re += a[i];
                    }
                    return f;
                },
                else => {
                    @compileError("inconsistent polynomial and variable type");
                },
            }
        },

        Complex(f32), Complex(f64) => {
            f = F.init(0, 0);
            switch (X) {
                f32, f64 => {
                    var i: usize = 0;
                    while (i < a.len) : (i += 1) {
                        f = cmath.add(a[i], cmath.scale(f, x));
                    }
                    return f;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    while (i < a.len) : (i += 1) {
                        f = cmath.add(a[i], cmath.mul(f, x));
                    }
                    return f;
                },
                else => {
                    @compileError("inconsistent polynomial and variable type");
                },
            }
        },
        else => {
            @compileError("inconsistent polynomial and variable type");
        },
    }
}

//---------------------------------------------------------------

const eps = 1e-4;

test "\t math,\t add - real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, 3);
        a[0] = 1.1;
        a[1] = 2.2;
        a[2] = -3.1;

        const b = try allocator.alloc(R, 4);
        b[0] = 1.1;
        b[1] = 2.2;
        b[2] = 3.3;
        b[3] = 1.2;

        const c = try allocator.alloc(R, 4);

        try add(R, c, a, b);

        try std.testing.expectApproxEqAbs(@as(R, 2.2), c[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.4), c[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.2), c[2], eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.2), c[3], eps);
    }
}

test "\t math,\t add - complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(C, 3);
        a[0].re = 1.1;
        a[0].im = 10.1;
        a[1].re = 2.2;
        a[1].im = 20.2;
        a[2].re = -3.1;
        a[2].im = -30.1;

        const b = try allocator.alloc(C, 4);
        b[0].re = 1.1;
        b[0].im = 10.1;
        b[1].re = 2.2;
        b[1].im = 20.2;
        b[2].re = 3.3;
        b[2].im = 30.3;
        b[3].re = 1.2;
        b[3].im = 10.2;

        const c = try allocator.alloc(C, 4);

        try add(C, c, a, b);

        try std.testing.expectApproxEqAbs(@as(R, 2.2), c[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 20.2), c[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.4), c[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 40.4), c[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.2), c[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.2), c[2].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.2), c[3].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 10.2), c[3].im, eps);
    }
}

test "\t math,\t conv - real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, 3);
        a[0] = 1.1;
        a[1] = 2.2;
        a[2] = -3.3;

        const b = try allocator.alloc(R, 4);
        b[0] = 1.1;
        b[1] = 2.2;
        b[2] = 3.3;
        b[3] = -4.4;

        const c = try allocator.alloc(R, 6);

        try conv(R, c, a, b);

        try std.testing.expectApproxEqAbs(@as(R, 1.21), c[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.84), c[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.84), c[2], eps);
        try std.testing.expectApproxEqAbs(@as(R, -4.84), c[3], eps);
        try std.testing.expectApproxEqAbs(@as(R, -20.57), c[4], eps);
        try std.testing.expectApproxEqAbs(@as(R, 14.52), c[5], eps);
    }
}

test "\t math,\t conv - complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(C, 3);
        a[0].re = 1.1;
        a[0].im = 1.0;
        a[1].re = 2.2;
        a[1].im = 3.0;
        a[2].re = -3.3;
        a[2].im = 4.0;

        const b = try allocator.alloc(C, 4);
        b[0].re = 1.1;
        b[0].im = 2.0;
        b[1].re = 2.2;
        b[1].im = 1.0;
        b[2].re = 3.3;
        b[2].im = 3.0;
        b[3].re = -4.4;
        b[3].im = -3.0;

        const c = try allocator.alloc(C, 6);

        try conv(C, c, a, b);

        try std.testing.expectApproxEqAbs(@as(R, -0.79), c[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 3.30), c[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -2.16), c[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 11.00), c[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -9.16), c[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 13.20), c[2].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -14.84), c[3].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 14.30), c[3].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -23.57), c[4].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -16.50), c[4].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 26.52), c[5].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -7.70), c[5].im, eps);
    }
}

test "\t math,\t evalDeriv - real/real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try allocator.alloc(R, 4);
        a[0] = 1.1;
        a[1] = 2.2;
        a[2] = -3.3;
        a[3] = 1.2;

        const x: R = -4.4;
        const f = evalDeriv(R, R, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 100.935999999), f, eps);
    }
}

test "\t math,\t evalDeriv - real/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, 4);
        a[0] = 1.1;
        a[1] = 2.2;
        a[2] = -3.3;
        a[3] = 1.2;

        const x = C.init(-0.2, 0.1);
        const f = evalDeriv(R, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 3.6280), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.804), f.im, eps);
    }
}

test "\t math,\t evalDeriv - complex/real \n" {
    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(C, 4);
        a[0] = C.init(1.20, 2.31);
        a[1] = C.init(-0.11, 2.70);
        a[2] = C.init(1.20, -1.33);
        a[3] = C.init(-0.40, 2.34);

        const x: R = -4.4;
        const f = evalDeriv(C, R, a, x);

        try std.testing.expectApproxEqAbs(@as(R, -33.90200), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 150.31120), f.im, eps);
    }
}

test "\t math,\t evalDeriv - complex/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(C, 4);
        a[0] = C.init(1.20, 2.31);
        a[1] = C.init(-0.11, 2.70);
        a[2] = C.init(1.20, -1.33);
        a[3] = C.init(-0.40, 2.34);

        const x = C.init(-0.2, 0.1);
        const f = evalDeriv(C, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, -0.0791999999), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 3.7306000000), f.im, eps);
    }
}

test "\t math,\t evalDerivRev - real/real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, 4);
        a[0] = 1.1;
        a[1] = 2.2;
        a[2] = -3.3;
        a[3] = 1.2;

        const x: R = -4.4;
        const f = evalDerivRev(R, R, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 41.2280), f, eps);
    }
}

test "\t math,\t evalDerivRev - real/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, 4);
        a[0] = 1.1;
        a[1] = 2.2;
        a[2] = -3.3;
        a[3] = 1.2;

        const x = C.init(-0.2, 0.1);
        const f = evalDerivRev(R, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, -4.0810), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.308000), f.im, eps);
    }
}

test "\t math,\t evalDerivRev - complex/real \n" {
    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(C, 4);
        a[0] = C.init(1.20, 2.31);
        a[1] = C.init(-0.11, 2.70);
        a[2] = C.init(1.20, -1.33);
        a[3] = C.init(-0.40, 2.34);

        const x: R = -4.4;
        const f = evalDerivRev(C, R, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 71.864), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 109.0748), f.im, eps);
    }
}

test "\t math,\t evalDerivRev - complex/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(C, 4);
        a[0] = C.init(1.20, 2.31);
        a[1] = C.init(-0.11, 2.70);
        a[2] = C.init(1.20, -1.33);
        a[3] = C.init(-0.40, 2.34);

        const x = C.init(-0.2, 0.1);
        const f = evalDerivRev(C, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 1.08920), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -2.36810), f.im, eps);
    }
}

test "\t math,\t eval - real/real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, 4);
        a[0] = 1.1;
        a[1] = 2.2;
        a[2] = -3.3;
        a[3] = 1.2;

        const x: R = -4.4;
        const f = eval(R, R, a, x);

        try std.testing.expectApproxEqAbs(@as(R, -174.68880000), f, eps);
    }
}

test "\t math,\t eval - real/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, 4);
        a[0] = 1.1;
        a[1] = 2.2;
        a[2] = -3.3;
        a[3] = 1.2;

        const x = C.init(-4.4, 1.0);
        const f = eval(R, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, -155.54880), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 99.736), f.im, eps);
    }
}

test "\t math,\t eval - complex/real \n" {
    inline for (.{ f32, f64 }) |R| {
        // inline for (.{f64}) |R| {

        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(C, 4);
        a[0] = C.init(1.20, 2.31);
        a[1] = C.init(-0.11, 2.70);
        a[2] = C.init(1.20, -1.33);
        a[3] = C.init(-0.40, 2.34);

        const x: R = -4.4;
        const f = eval(C, R, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 58.9896), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -234.6493600), f.im, eps);
    }
}

test "\t math,\t eval - complex/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(C, 4);
        a[0] = C.init(1.20, 2.31);
        a[1] = C.init(-0.11, 2.70);
        a[2] = C.init(1.20, -1.33);
        a[3] = C.init(-0.40, 2.34);

        const x = C.init(-0.2, 0.1);
        const f = eval(C, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 0.9098599999), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.6620200000), f.im, eps);
    }
}

test "\t math,\t evalRev - real/real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, 4);
        a[0] = 1.1;
        a[1] = 2.2;
        a[2] = -3.3;
        a[3] = 1.2;

        const x: R = -4.4;
        const f = evalRev(R, R, a, x);

        try std.testing.expectApproxEqAbs(f, -35.3904000, eps);
    }
}

test "\t math,\t evalRev - real/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, 4);
        a[0] = 1.1;
        a[1] = 2.2;
        a[2] = -3.3;
        a[3] = 1.2;

        const x = C.init(-0.2, 0.1);
        const f = evalRev(R, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 1.9238), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.40590), f.im, eps);
    }
}

test "\t math,\t evalRev - complex/real \n" {
    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(C, 4);
        a[0] = C.init(1.20, 2.31);
        a[1] = C.init(-0.11, 2.70);
        a[2] = C.init(1.20, -1.33);
        a[3] = C.init(-0.40, 2.34);

        const x: R = -4.4;
        const f = evalRev(C, R, a, x);

        try std.testing.expectApproxEqAbs(@as(R, -110.0304), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -136.31104), f.im, eps);
    }
}

test "\t math,\t evalRev - complex/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(C, 4);
        a[0] = C.init(1.20, 2.31);
        a[1] = C.init(-0.11, 2.70);
        a[2] = C.init(1.20, -1.33);
        a[3] = C.init(-0.40, 2.34);

        const x = C.init(-0.2, 0.1);
        const f = evalRev(C, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, -0.43011), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 2.819980), f.im, eps);
    }
}
