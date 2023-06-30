const std = @import("std");
const print = std.debug.print;
const math = std.math;
const Complex = std.math.complex.Complex;

const zpoly = @import("../zpoly.zig");
const Polynomial = zpoly.Polynomial;

const EitherCmpx = zpoly.EitherCmpx;
const ValueType = zpoly.ValueType;

pub fn add(comptime T: type, a: *Polynomial(T), f: Polynomial(T), g: Polynomial(T)) !void {
    //
    // add polynomials f(x) and g(x) of potentially different orders.
    //
    //    a(x) = f(x) + g(x)
    //
    //    if order of f(x) is p (has p+1 coef), and order of g(x) is
    //    q (has q+1 coef), then a(x) will have order max(p, q)
    //    (has max(p+1, q+1) coeff.
    //
    //
    // output polynomial a(x) must have sufficient space to store result.

    var min_len: usize = @min(f.val.len, g.val.len);
    var max_len: usize = @max(f.val.len, g.val.len);

    if (a.maxlen < max_len) {
        return error.Insufficient_Len_to_Store_Result;
    }
    a.val.len = max_len;

    switch (T) {
        f32, f64 => {
            var i: usize = 0;
            while (i < min_len) : (i += 1) {
                a.val[i] = f.val[i] + g.val[i];
            }

            if (f.val.len > g.val.len) {
                while (i < f.val.len) : (i += 1) {
                    a.val[i] = f.val[i];
                }
            } else if (g.val.len > f.val.len) {
                while (i < g.val.len) : (i += 1) {
                    a.val[i] = g.val[i];
                }
            }
        },
        Complex(f32), Complex(f64) => {
            var i: usize = 0;
            while (i < min_len) : (i += 1) {
                a.val[i].re = f.val[i].re + g.val[i].re;
                a.val[i].im = f.val[i].im + g.val[i].im;
            }

            if (f.val.len > g.val.len) {
                while (i < f.val.len) : (i += 1) {
                    a.val[i].re = f.val[i].re;
                    a.val[i].im = f.val[i].im;
                }
            } else if (g.val.len > f.val.len) {
                while (i < g.val.len) : (i += 1) {
                    a.val[i].re = g.val[i].re;
                    a.val[i].im = g.val[i].im;
                }
            }
        },
        else => {
            @compileError("unknown polynomial type");
        },
    }
}

pub fn evalDeriv(comptime A: type, comptime X: type, a: Polynomial(A), x: X) EitherCmpx(A, X) {
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
                    while (i < a.val.len) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.val.len - i)));
                        g = n_tmp * a.val[a.val.len - i] + x * g;
                    }
                    return g;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 1;
                    var g_tmp: G = undefined;
                    g = G.init(0, 0);
                    while (i < a.val.len) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.val.len - i)));
                        var a_tmp: A = a.val[a.val.len - i];
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
                    var i: usize = 1;
                    var g_tmp: G = undefined;
                    var a_tmp: G = undefined;
                    g = G.init(0, 0);
                    while (i < a.val.len) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.val.len - i)));
                        a_tmp.re = a.val[a.val.len - i].re;
                        a_tmp.im = a.val[a.val.len - i].im;

                        g_tmp.re = g.re;
                        g_tmp.im = g.im;

                        g.re = n_tmp * a_tmp.re + (x * g_tmp.re);
                        g.im = n_tmp * a_tmp.im + (x * g_tmp.im);
                    }
                    return g;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 1;
                    var g_tmp: G = undefined;
                    var a_tmp: G = undefined;
                    g = G.init(0, 0);
                    while (i < a.val.len) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.val.len - i)));
                        a_tmp.re = a.val[a.val.len - i].re;
                        a_tmp.im = a.val[a.val.len - i].im;

                        g_tmp.re = g.re;
                        g_tmp.im = g.im;

                        g.re = n_tmp * a_tmp.re + (x.re * g_tmp.re - x.im * g_tmp.im);
                        g.im = n_tmp * a_tmp.im + (x.im * g_tmp.re + x.re * g_tmp.im);
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

pub fn evalDerivRev(comptime A: type, comptime X: type, a: Polynomial(A), x: X) EitherCmpx(A, X) {
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
                    while (i < a.val.len - 1) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.val.len - i - 1)));
                        g = n_tmp * a.val[i] + x * g;
                    }
                    return g;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    var g_tmp: G = undefined;
                    g = G.init(0, 0);
                    while (i < a.val.len - 1) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.val.len - i - 1)));
                        var a_tmp: A = a.val[i];
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
                    var g_tmp: G = undefined;
                    var a_tmp: G = undefined;
                    g = G.init(0, 0);
                    while (i < a.val.len - 1) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.val.len - i - 1)));
                        a_tmp.re = a.val[i].re;
                        a_tmp.im = a.val[i].im;

                        g_tmp.re = g.re;
                        g_tmp.im = g.im;

                        g.re = n_tmp * a_tmp.re + (x * g_tmp.re);
                        g.im = n_tmp * a_tmp.im + (x * g_tmp.im);
                    }
                    return g;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    var g_tmp: G = undefined;
                    var a_tmp: G = undefined;
                    g = G.init(0, 0);
                    while (i < a.val.len - 1) : (i += 1) {
                        n_tmp = @as(R, @floatFromInt((a.val.len - i - 1)));
                        a_tmp.re = a.val[i].re;
                        a_tmp.im = a.val[i].im;

                        g_tmp.re = g.re;
                        g_tmp.im = g.im;

                        g.re = n_tmp * a_tmp.re + (x.re * g_tmp.re - x.im * g_tmp.im);
                        g.im = n_tmp * a_tmp.im + (x.im * g_tmp.re + x.re * g_tmp.im);
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

pub fn eval(comptime A: type, comptime X: type, a: Polynomial(A), x: X) EitherCmpx(A, X) {
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
                    while (i < a.val.len) : (i += 1) {
                        f = a.val[a.val.len - 1 - i] + x * f;
                    }
                    return f;
                },
                Complex(f32), Complex(f64) => {
                    var a_tmp: A = undefined;
                    var f_tmp: F = undefined;
                    f.re = 0;
                    f.im = 0;
                    var i: usize = 0;
                    while (i < a.val.len) : (i += 1) {
                        a_tmp = a.val[a.val.len - 1 - i];
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
            var a_tmp: A = undefined;
            var f_tmp: F = undefined;
            f.re = 0;
            f.im = 0;
            switch (X) {
                f32, f64 => {
                    var i: usize = 0;
                    while (i < a.val.len) : (i += 1) {
                        a_tmp.re = a.val[a.val.len - 1 - i].re;
                        a_tmp.im = a.val[a.val.len - 1 - i].im;
                        f_tmp.re = f.re;
                        f_tmp.im = f.im;
                        f.re = a_tmp.re + (x * f_tmp.re);
                        f.im = a_tmp.im + (x * f_tmp.im);
                    }
                    return f;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    while (i < a.val.len) : (i += 1) {
                        a_tmp.re = a.val[a.val.len - 1 - i].re;
                        a_tmp.im = a.val[a.val.len - 1 - i].im;
                        f_tmp.re = f.re;
                        f_tmp.im = f.im;
                        f.re = a_tmp.re + (x.re * f_tmp.re - x.im * f_tmp.im);
                        f.im = a_tmp.im + (x.im * f_tmp.re + x.re * f_tmp.im);
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

pub fn mult(comptime T: type, a: *Polynomial(T), f: Polynomial(T), g: Polynomial(T)) !void {
    // multiply two polynomials
    //
    //            a(x) = f(x) g(x)
    //
    // a(x) is of order p+q (having p+q+1 coefs) where order of f(x)=p,
    // (having p+1 coefs) and order of g(x) = q (having q+1 coefs).
    //
    // a must have enough space to store the result

    if (a.maxlen < f.val.len + g.val.len + 1) {
        return error.Insufficient_Space_to_Store_Result;
    }
    a.val.len = f.val.len + g.val.len - 1;

    var n: usize = 0;

    switch (T) {
        f32, f64 => {
            while (n < a.val.len) : (n += 1) {
                a.val[n] = 0;
                var k: usize = 0;
                while (k < f.val.len) : (k += 1) {
                    if ((n >= k) and (n - k < g.val.len)) {
                        a.val[n] += f.val[k] * g.val[n - k];
                    }
                }
            }
        },
        Complex(f32), Complex(f64) => {
            while (n < a.val.len) : (n += 1) {
                a.val[n].re = 0;
                a.val[n].im = 0;

                var k: usize = 0;
                while (k < f.val.len) : (k += 1) {
                    if ((n >= k) and (n - k < g.val.len)) {
                        a.val[n].re += f.val[k].re * g.val[n - k].re - f.val[k].im * g.val[n - k].im;
                        a.val[n].im += f.val[k].re * g.val[n - k].im + f.val[k].im * g.val[n - k].re;
                    }
                }
            }
        },

        else => {
            @compileError("unknown polynomial type");
        },
    }
}

pub fn evalRev(comptime A: type, comptime X: type, a: Polynomial(A), x: X) EitherCmpx(A, X) {
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
                    while (i < a.val.len) : (i += 1) {
                        f = a.val[i] + x * f;
                    }
                    return f;
                },
                Complex(f32), Complex(f64) => {
                    var a_tmp: A = undefined;
                    var f_tmp: F = undefined;
                    f.re = 0;
                    f.im = 0;
                    var i: usize = 0;
                    while (i < a.val.len) : (i += 1) {
                        a_tmp = a.val[i];
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
            var a_tmp: A = undefined;
            var f_tmp: F = undefined;
            f.re = 0;
            f.im = 0;
            switch (X) {
                f32, f64 => {
                    var i: usize = 0;
                    while (i < a.val.len) : (i += 1) {
                        a_tmp.re = a.val[i].re;
                        a_tmp.im = a.val[i].im;
                        f_tmp.re = f.re;
                        f_tmp.im = f.im;
                        f.re = a_tmp.re + (x * f_tmp.re);
                        f.im = a_tmp.im + (x * f_tmp.im);
                    }
                    return f;
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    while (i < a.val.len) : (i += 1) {
                        a_tmp.re = a.val[i].re;
                        a_tmp.im = a.val[i].im;
                        f_tmp.re = f.re;
                        f_tmp.im = f.im;
                        f.re = a_tmp.re + (x.re * f_tmp.re - x.im * f_tmp.im);
                        f.im = a_tmp.im + (x.im * f_tmp.re + x.re * f_tmp.im);
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

        var a = try Polynomial(R).init(allocator, 4, 3);
        a.val[0] = 1.1;
        a.val[1] = 2.2;
        a.val[2] = -3.1;

        var b = try Polynomial(R).init(allocator, 6, 4);
        b.val[0] = 1.1;
        b.val[1] = 2.2;
        b.val[2] = 3.3;
        b.val[3] = 1.2;

        var c = try Polynomial(R).init(allocator, 10, 10);

        try add(R, &c, a, b);

        try std.testing.expectEqual(@as(usize, 4), c.val.len);
        try std.testing.expectApproxEqAbs(@as(R, 2.2), c.val[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.4), c.val[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.2), c.val[2], eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.2), c.val[3], eps);
    }
}

test "\t math,\t add - complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try Polynomial(C).init(allocator, 4, 3);
        a.val[0].re = 1.1;
        a.val[0].im = 10.1;
        a.val[1].re = 2.2;
        a.val[1].im = 20.2;
        a.val[2].re = -3.1;
        a.val[2].im = -30.1;

        var b = try Polynomial(C).init(allocator, 6, 4);
        b.val[0].re = 1.1;
        b.val[0].im = 10.1;
        b.val[1].re = 2.2;
        b.val[1].im = 20.2;
        b.val[2].re = 3.3;
        b.val[2].im = 30.3;
        b.val[3].re = 1.2;
        b.val[3].im = 10.2;

        var c = try Polynomial(C).init(allocator, 10, 10);

        try add(C, &c, a, b);

        try std.testing.expectEqual(c.val.len, 4);
        try std.testing.expectApproxEqAbs(@as(R, 2.2), c.val[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 20.2), c.val[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.4), c.val[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 40.4), c.val[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.2), c.val[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.2), c.val[2].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.2), c.val[3].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 10.2), c.val[3].im, eps);
    }
}

test "\t math,\t evalDeriv - real/real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try Polynomial(R).init(allocator, 4, 4);
        a.val[0] = 1.1;
        a.val[1] = 2.2;
        a.val[2] = -3.3;
        a.val[3] = 1.2;

        var x: R = -4.4;
        var f = evalDeriv(R, R, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 100.935999999), f, eps);
    }
}

test "\t math,\t evalDeriv - real/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try Polynomial(R).init(allocator, 4, 4);
        a.val[0] = 1.1;
        a.val[1] = 2.2;
        a.val[2] = -3.3;
        a.val[3] = 1.2;

        var x = C.init(-0.2, 0.1);
        var f = evalDeriv(R, C, a, x);

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

        var a = try Polynomial(C).init(allocator, 4, 4);
        a.val[0] = C.init(1.20, 2.31);
        a.val[1] = C.init(-0.11, 2.70);
        a.val[2] = C.init(1.20, -1.33);
        a.val[3] = C.init(-0.40, 2.34);

        var x: R = -4.4;
        var f = evalDeriv(C, R, a, x);

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

        var a = try Polynomial(C).init(allocator, 4, 4);
        a.val[0] = C.init(1.20, 2.31);
        a.val[1] = C.init(-0.11, 2.70);
        a.val[2] = C.init(1.20, -1.33);
        a.val[3] = C.init(-0.40, 2.34);

        var x = C.init(-0.2, 0.1);
        var f = evalDeriv(C, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, -0.0791999999), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 3.7306000000), f.im, eps);
    }
}

test "\t math,\t evalDerivRev - real/real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try Polynomial(R).init(allocator, 4, 4);
        a.val[0] = 1.1;
        a.val[1] = 2.2;
        a.val[2] = -3.3;
        a.val[3] = 1.2;

        var x: R = -4.4;
        var f = evalDerivRev(R, R, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 41.2280), f, eps);
    }
}

test "\t math,\t evalDerivRev - real/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try Polynomial(R).init(allocator, 4, 4);
        a.val[0] = 1.1;
        a.val[1] = 2.2;
        a.val[2] = -3.3;
        a.val[3] = 1.2;

        var x = C.init(-0.2, 0.1);
        var f = evalDerivRev(R, C, a, x);

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

        var a = try Polynomial(C).init(allocator, 4, 4);
        a.val[0] = C.init(1.20, 2.31);
        a.val[1] = C.init(-0.11, 2.70);
        a.val[2] = C.init(1.20, -1.33);
        a.val[3] = C.init(-0.40, 2.34);

        var x: R = -4.4;
        var f = evalDerivRev(C, R, a, x);

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

        var a = try Polynomial(C).init(allocator, 4, 4);
        a.val[0] = C.init(1.20, 2.31);
        a.val[1] = C.init(-0.11, 2.70);
        a.val[2] = C.init(1.20, -1.33);
        a.val[3] = C.init(-0.40, 2.34);

        var x = C.init(-0.2, 0.1);
        var f = evalDerivRev(C, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 1.08920), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -2.36810), f.im, eps);
    }
}

test "\t math,\t eval - real/real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try Polynomial(R).init(allocator, 4, 4);
        a.val[0] = 1.1;
        a.val[1] = 2.2;
        a.val[2] = -3.3;
        a.val[3] = 1.2;

        var x: R = -4.4;
        var f = eval(R, R, a, x);

        try std.testing.expectApproxEqAbs(@as(R, -174.68880000), f, eps);
    }
}

test "\t math,\t eval - real/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try Polynomial(R).init(allocator, 4, 4);
        a.val[0] = 1.1;
        a.val[1] = 2.2;
        a.val[2] = -3.3;
        a.val[3] = 1.2;

        var x = C.init(-4.4, 1.0);
        var f = eval(R, C, a, x);

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

        var a = try Polynomial(C).init(allocator, 4, 4);
        a.val[0] = C.init(1.20, 2.31);
        a.val[1] = C.init(-0.11, 2.70);
        a.val[2] = C.init(1.20, -1.33);
        a.val[3] = C.init(-0.40, 2.34);

        var x: R = -4.4;
        var f = eval(C, R, a, x);

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

        var a = try Polynomial(C).init(allocator, 4, 4);
        a.val[0] = C.init(1.20, 2.31);
        a.val[1] = C.init(-0.11, 2.70);
        a.val[2] = C.init(1.20, -1.33);
        a.val[3] = C.init(-0.40, 2.34);

        var x = C.init(-0.2, 0.1);
        var f = eval(C, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, 0.9098599999), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.6620200000), f.im, eps);
    }
}

test "\t math,\t mult - real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try Polynomial(R).init(allocator, 4, 3);
        a.val[0] = 1.1;
        a.val[1] = 2.2;
        a.val[2] = -3.3;

        var b = try Polynomial(R).init(allocator, 6, 4);
        b.val[0] = 1.1;
        b.val[1] = 2.2;
        b.val[2] = 3.3;
        b.val[3] = -4.4;

        var c = try Polynomial(R).init(allocator, 10, 10);

        try mult(R, &c, a, b);

        try std.testing.expectEqual(c.val.len, 6);
        try std.testing.expectApproxEqAbs(@as(R, 1.21), c.val[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.84), c.val[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.84), c.val[2], eps);
        try std.testing.expectApproxEqAbs(@as(R, -4.84), c.val[3], eps);
        try std.testing.expectApproxEqAbs(@as(R, -20.57), c.val[4], eps);
        try std.testing.expectApproxEqAbs(@as(R, 14.52), c.val[5], eps);
    }
}

test "\t math,\t mult - complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try Polynomial(C).init(allocator, 4, 3);
        a.val[0].re = 1.1;
        a.val[0].im = 1.0;
        a.val[1].re = 2.2;
        a.val[1].im = 3.0;
        a.val[2].re = -3.3;
        a.val[2].im = 4.0;

        var b = try Polynomial(C).init(allocator, 6, 4);
        b.val[0].re = 1.1;
        b.val[0].im = 2.0;
        b.val[1].re = 2.2;
        b.val[1].im = 1.0;
        b.val[2].re = 3.3;
        b.val[2].im = 3.0;
        b.val[3].re = -4.4;
        b.val[3].im = -3.0;

        var c = try Polynomial(C).init(allocator, 10, 10);

        try mult(C, &c, a, b);

        try std.testing.expectEqual(c.val.len, 6);
        try std.testing.expectApproxEqAbs(@as(R, -0.79), c.val[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 3.30), c.val[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -2.16), c.val[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 11.00), c.val[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -9.16), c.val[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 13.20), c.val[2].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -14.84), c.val[3].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 14.30), c.val[3].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -23.57), c.val[4].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -16.50), c.val[4].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 26.52), c.val[5].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -7.70), c.val[5].im, eps);
    }
}

test "\t math,\t evalRev - real/real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try Polynomial(R).init(allocator, 4, 4);
        a.val[0] = 1.1;
        a.val[1] = 2.2;
        a.val[2] = -3.3;
        a.val[3] = 1.2;

        var x: R = -4.4;
        var f = evalRev(R, R, a, x);

        try std.testing.expectApproxEqAbs(f, -35.3904000, eps);
    }
}

test "\t math,\t evalRev - real/complex \n" {
    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try Polynomial(R).init(allocator, 4, 4);
        a.val[0] = 1.1;
        a.val[1] = 2.2;
        a.val[2] = -3.3;
        a.val[3] = 1.2;

        var x = C.init(-0.2, 0.1);
        var f = evalRev(R, C, a, x);

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

        var a = try Polynomial(C).init(allocator, 4, 4);
        a.val[0] = C.init(1.20, 2.31);
        a.val[1] = C.init(-0.11, 2.70);
        a.val[2] = C.init(1.20, -1.33);
        a.val[3] = C.init(-0.40, 2.34);

        var x: R = -4.4;
        var f = evalRev(C, R, a, x);

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

        var a = try Polynomial(C).init(allocator, 4, 4);
        a.val[0] = C.init(1.20, 2.31);
        a.val[1] = C.init(-0.11, 2.70);
        a.val[2] = C.init(1.20, -1.33);
        a.val[3] = C.init(-0.40, 2.34);

        var x = C.init(-0.2, 0.1);
        var f = evalRev(C, C, a, x);

        try std.testing.expectApproxEqAbs(@as(R, -0.43011), f.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 2.819980), f.im, eps);
    }
}
