const std = @import("std");
const Allocator = std.mem.Allocator;

const Child = std.meta.Child;
const bytesAsSlice = std.mem.bytesAsSlice;

const print = std.debug.print;
const math = std.math;
const Complex = std.math.Complex;

const zpoly = @import("../zpoly.zig");
const cmath = zpoly.cmath;

const fabs = zpoly.fabs;
const eval = zpoly.eval;
const evalDeriv = zpoly.evalDeriv;
const evalRev = zpoly.evalRev;
const evalDerivRev = zpoly.evalDerivRev;
const ValueType = zpoly.ValueType;
const ToCmpx = zpoly.ToCmpx;

fn allocStructSlices(allocator: Allocator, X: anytype, args: anytype) ![]align(@alignOf(usize)) u8 {

    // detemine n_bytes needed (add extra slop for alignments)
    var n_bytes: usize = 0;
    inline for (std.meta.fields(@TypeOf(X.*)), 0..) |f, i| {
        n_bytes += @sizeOf(f.type) + args[i] * @sizeOf(Child(f.type));
        n_bytes += 4 * @sizeOf(usize);
    }

    // allocate the raw bytes
    var raw_bytes: []align(@alignOf(usize)) u8 =
        try allocator.allocWithOptions(u8, n_bytes, @alignOf(usize), null);

    var start: usize = 0;
    var end: usize = undefined;
    var f_bytes: usize = 0;

    inline for (std.meta.fields(@TypeOf(X.*)), 0..) |f, i| {

        // calculate end for this struct field
        const child_type: type = Child(f.type);
        f_bytes = @sizeOf(f.type) + args[i] * @sizeOf(child_type);
        end = start + f_bytes;

        // partition raw into slice for field, and set length
        @field(X.*, @typeInfo(@TypeOf(X.*)).Struct.fields[i].name) =
            @alignCast(bytesAsSlice(child_type, raw_bytes[start..end]));

        @field(X.*, @typeInfo(@TypeOf(X.*)).Struct.fields[i].name).len = args[i];

        // calculate start for nextnext struct field
        start = end + @sizeOf(usize) - @rem(end, @sizeOf(usize));
    }

    // return original raw_bytes allocated to allow free later
    return raw_bytes;
}

// This code follows D.A. Bini, "Numerical computation of polynomial zeros by means of
// Aberth's method", Numerical Algorithms, Feb 1996.
//
// The notation and variable names used here agree with article wherever possible.

fn stop_criterion(comptime P: type, comptime C: type, a: []P, s_tilde: []ValueType(P), x: C) bool {

    // Bini equation 15

    const R: type = ValueType(P);
    var mu: R = undefined;

    switch (R) {
        f32 => {
            mu = 1e-6;
        },
        f64 => {
            mu = 1e-14;
        },
        else => {
            @compileError("unexpected type");
        },
    }

    const lhs: R = fabs(C, eval(P, C, a, x));
    const rhs: R = mu * eval(R, R, s_tilde, fabs(C, x)); // this is strictly non-negative

    if (lhs < rhs) {
        return true;
    } else {
        return false;
    }
}

fn product(comptime T: type, x: T, y: T) T {
    switch (T) {
        f32, f64 => {
            return x * y;
        },
        Complex(f32), Complex(f64) => {
            return cmath.mul(x, y);
        },
        else => {
            @compileError("product for data type not implement");
        },
    }
}

fn quotient(comptime T: type, num: T, den: T) T {
    switch (T) {
        f32, f64 => {
            return num / den;
        },
        Complex(f32), Complex(f64) => {
            const tmp: ValueType(T) = 1 / (den.re * den.re + den.im * den.im);

            var quo: T = undefined;
            quo.re = tmp * (num.re * den.re + num.im * den.im);
            quo.im = tmp * (-num.re * den.im + num.im * den.re);
            return quo;
        },
        else => {
            @compileError("reciprocal for data type not implement");
        },
    }
}

fn recip(comptime T: type, x: T) T {
    switch (T) {
        f32, f64 => {
            return 1 / x;
        },
        Complex(f32), Complex(f64) => {
            return cmath.recip(x);
        },
        else => {
            @compileError("reciprocal for data type not implement");
        },
    }
}

pub fn roots(comptime A: type, allocator: Allocator, a: []A, x: []ToCmpx(A)) !void {
    var idx_high: usize = a.len;
    while (idx_high != 0) {
        idx_high -= 1;
        if (fabs(A, a[idx_high]) != 0) {
            break;
        }
    }

    var idx_low: usize = 0;
    while (idx_low != idx_high) : (idx_low += 1) {
        if (fabs(A, a[idx_low]) != 0) {
            break;
        }
    }

    var i: usize = 0;
    while (i < x.len) : (i += 1) {
        x[i].re = 0;
        x[i].im = 0;
    }

    try roots_internal(A, allocator, a[idx_low .. idx_high + 1], x[0 .. idx_high - idx_low]);
}

fn roots_internal(comptime A: type, allocator: Allocator, a: []A, x: []ToCmpx(A)) !void {
    const R = ValueType(A);

    // root finding algorithm requires non-zero leading and trailing coefficients,
    // (e.g. a_n * a_0 != 0).
    if (fabs(A, a[a.len - 1]) == 0) {
        return error.LeadingCoefficientIsZero;
    } else if (fabs(A, a[0]) == 0) {
        return error.TrailingCoefficientIsZero;
    }

    var s_tilde: []R = undefined;

    const Workspace = struct {
        k: []usize,
        hull: []R,
        u: []R,
        s_tilde_tmp: []R,
        stop_flags: []bool,
    };

    const len: usize = a.len;
    var workspace: Workspace = undefined;
    const raw_bytes = try allocStructSlices(allocator, &workspace, .{ len, len, len, len, len });

    // see lines below equation 14 in Bini
    // s_tilde has strictly real coefficients, evaluated on strictly real arguments

    for (a, 0..) |aval, i| {
        workspace.s_tilde_tmp[i] = fabs(A, aval) * @as(R, @floatFromInt((4 * i + 1)));
    }
    s_tilde = workspace.s_tilde_tmp;

    for (workspace.stop_flags, 0..) |_, i| {
        workspace.stop_flags[i] = false;
    }

    upperConvexHull(A, a, &workspace.k, &workspace.hull);
    computeModulii(A, a, workspace.k, workspace.u);
    initialRootEst(R, workspace.k, workspace.u, x);

    const max_iter = 20;
    aberth(A, max_iter, a, s_tilde, workspace.stop_flags, x);
    allocator.free(raw_bytes);
}

pub fn upperConvexHull(comptime A: type, a: []A, k: *[]usize, hull: *[]ValueType(A)) void {
    const R = ValueType(A);

    const tiny_neg = -1e-20;

    k.ptr[0] = 0;
    k.ptr[1] = 1;

    if (fabs(A, a[0]) == 0) {
        hull.ptr[0] = tiny_neg;
    } else {
        hull.ptr[0] = math.log(R, 2.0, fabs(A, a[0]));
    }

    if (fabs(A, a[1]) == 0) {
        hull.ptr[1] = tiny_neg;
    } else {
        hull.ptr[1] = math.log(R, 2.0, fabs(A, a[1]));
    }

    var len: usize = 2;
    var i: usize = 2;
    while (i < a.len) : (i += 1) {

        // append to lists
        len += 1;
        k.ptr[len - 1] = i;

        if (fabs(A, a[i]) == 0) {
            hull.ptr[len - 1] = tiny_neg;
        } else {
            hull.ptr[len - 1] = math.log(R, 2.0, fabs(A, a[i]));
        }

        while (len > 2) {

            // to detect a "right turn" use sign of vector cross-product,
            // visualize three points on the plane,  a, b, c, as two vectors:
            // 1st vector joins a w/ b (called "ab"), has components ab_x, and ab_y
            // 2nd vector joins b w/ c (called "bc"), etc.

            const ab_x: R = @as(R, @floatFromInt(k.ptr[len - 2] - k.ptr[len - 3]));
            const bc_x: R = @as(R, @floatFromInt(k.ptr[len - 1] - k.ptr[len - 2]));

            const ab_y: R = hull.ptr[len - 2] - hull.ptr[len - 3];
            const bc_y: R = hull.ptr[len - 1] - hull.ptr[len - 2];

            const cross: R = ab_x * bc_y - ab_y * bc_x;
            if (cross >= 0) {
                // if does not make a "right turn" then
                // delete middle point (point b) from list
                k.ptr[len - 2] = k.ptr[len - 1];
                hull.ptr[len - 2] = hull.ptr[len - 1];
                len -= 1;
            } else {
                break;
            }
        }
    }

    k.len = len;
    hull.len = len;
}

pub fn computeModulii(comptime A: type, a: []A, k: []usize, u: []ValueType(A)) void {
    const R = ValueType(A);

    const q: usize = k.len;
    var base: R = undefined;
    var exp: R = undefined;

    var i: usize = 1;
    while (i < q) : (i += 1) {
        base = fabs(A, a[k[i - 1]]) / fabs(A, a[k[i]]);
        exp = (1.0 / @as(R, @floatFromInt((k[i] - k[i - 1]))));
        u[i - 1] = math.pow(R, base, exp);
    }
    // NOTICE that (unlike in Bini's paper, which does not use index 0,
    // here we choose to use index 0 of u is used!
}

pub fn initialRootEst(comptime T: type, k: []usize, u: []T, r: []Complex(T)) void {
    const q: usize = k.len;
    const n: usize = k[k.len - 1];

    var alpha: T = undefined;
    var beta: T = undefined;
    const sigma = 0.7;
    var phi: T = undefined;

    var i: usize = 0;
    while (i < q - 1) : (i += 1) {
        var j: usize = 0;
        while (j < k[i + 1] - k[i]) : (j += 1) {
            alpha = 2 * math.pi * @as(T, @floatFromInt(j + 1)) / @as(T, @floatFromInt(k[i + 1] - k[i]));
            beta = (2 * math.pi * @as(T, @floatFromInt((i + 1)))) / @as(T, @floatFromInt(n));

            phi = alpha + beta + sigma;
            r[k[i] + j].re = u[i] * @cos(phi);
            r[k[i] + j].im = u[i] * @sin(phi);
        }
    }
}

pub fn aberth(comptime P: type, max_iter: usize, p: []P, s_tilde: []ValueType(P), stop_flags: []bool, x: []ToCmpx(P)) void {
    const C = ToCmpx(P);
    const R = ValueType(P);

    const n: usize = p.len - 1; // number of roots

    var px_dpx: C = undefined; // p(x) / p'(x)

    var update_den: C = undefined; // Aberth update denominator
    var update: C = undefined; // Aberth update

    var k: usize = 0;

    while (k < max_iter) : (k += 1) {
        var i: usize = 0;
        while (i < n) : (i += 1) {
            const x_mag: R = fabs(C, x[i]);

            if (x_mag <= 1.0) {
                px_dpx = quotient(C, eval(P, C, p, x[i]), evalDeriv(P, C, p, x[i]));
            } else {
                const gamma: C = recip(C, x[i]);
                const gamma_sqr: C = product(C, gamma, gamma);

                // tmp is p'_rev(gamma) / p_rev(gamma)    (eqn 13 in Bini paper)
                const tmp = quotient(C, evalDerivRev(P, C, p, gamma), evalRev(P, C, p, gamma));

                var term1: C = undefined;
                var term2: C = undefined;

                term1 = cmath.scale(gamma, @as(R, @floatFromInt(n)));

                term2 = product(C, gamma_sqr, tmp);

                px_dpx = cmath.sub(term1, term2);

                px_dpx = recip(C, px_dpx);
            }

            var sum_recip_diff_x = C.init(0, 0);
            var j: usize = 0;

            while (j < n) : (j += 1) {
                if (j != i) {
                    var diff_x = cmath.sub(x[i], x[j]);
                    diff_x = recip(C, diff_x);
                    sum_recip_diff_x = cmath.add(sum_recip_diff_x, diff_x);
                }
            }

            update_den = product(C, px_dpx, sum_recip_diff_x);
            update_den = cmath.sub(C.init(1, 0), update_den);

            update = quotient(C, px_dpx, update_den);

            // stopping criterion
            if (!stop_flags[i]) {
                stop_flags[i] = stop_criterion(P, C, p, s_tilde, x[i]);

                x[i] = cmath.sub(x[i], update);
            }
        }
    }
}

//-------------- TESTS BELOW -------------------------

test "\t roots,\t private fn \t allocStuctSlices\n" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const TestStruct = struct {
        a: []bool,
        b: []f64,
        c: []i32,
    };

    var ts: TestStruct = undefined;
    var raw_bytes: []align(@alignOf(usize)) u8 = undefined;
    raw_bytes = try allocStructSlices(allocator, &ts, .{ 5, 6, 7 });

    // print("@TypeOf(raw_bytes): {any}, raw_bytes.len: {d}\n", .{@TypeOf(raw_bytes), raw_bytes.len});
    // print("@TypeOf(ts.a): {any}, ts.a.len: {d}\n", .{@TypeOf(ts.a), ts.a.len});
    // print("@TypeOf(ts.b): {any}, ts.b.len: {d}\n", .{@TypeOf(ts.b), ts.b.len});
    // print("@TypeOf(ts.c): {any}, ts.c.len: {d}\n", .{@TypeOf(ts.c), ts.c.len});

    try std.testing.expectEqual(@as(usize, 225), raw_bytes.len);
    try std.testing.expectEqual(@as(usize, 5), ts.a.len);
    try std.testing.expectEqual(@as(usize, 6), ts.b.len);
    try std.testing.expectEqual(@as(usize, 7), ts.c.len);
}

test "\t roots,\t private fn \t product, real\n" {
    const eps = 1e-5;

    inline for (.{ f32, f64 }) |R| {
        const x: R = 10.0;
        const y: R = 2.2;
        var p: R = undefined;

        p = product(R, x, y);

        try std.testing.expectApproxEqAbs(@as(R, 22.0), p, eps);
    }
}

test "\t roots,\t private fn \t product, complex\n" {
    const eps = 1e-5;

    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        const x = C.init(0.2, -0.13);
        const y = C.init(1.5, 0.11);
        var p = C.init(undefined, undefined);

        p = product(C, x, y);

        try std.testing.expectApproxEqAbs(@as(R, 0.3143), p.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.17300), p.im, eps);
    }
}

test "\t roots,\t private fn \t quotient, real\n" {
    const eps = 1e-5;

    inline for (.{ f32, f64 }) |R| {
        const x: R = 10.0;
        const y: R = 2.2;
        var q: R = undefined;

        q = quotient(R, x, y);

        try std.testing.expectApproxEqAbs(@as(R, 4.5454545454), q, eps);
    }
}

test "\t roots,\t private fn \t quotient, complex\n" {
    const eps = 1e-5;

    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        const x = C.init(0.2, -0.13);
        const y = C.init(1.5, 0.11);
        var q = C.init(undefined, undefined);

        q = quotient(C, x, y);

        try std.testing.expectApproxEqAbs(@as(R, 0.1262985721), q.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.0959285561), q.im, eps);
    }
}

test "\t roots,\t private fn \t recip, real\n" {
    const eps = 1e-5;

    inline for (.{ f32, f64 }) |R| {
        const x: R = 10.0;
        var y: R = undefined;

        y = recip(R, x);

        try std.testing.expectApproxEqAbs(@as(R, 0.1), y, eps);
    }
}

test "\t roots,\t private fn \t recip, complex\n" {
    const eps = 1e-5;

    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);

        const x = C.init(0.20, -0.13);
        var y = C.init(undefined, undefined);

        y = recip(C, x);

        try std.testing.expectApproxEqAbs(@as(R, 3.514938488576), y.re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 2.284710017575), y.im, eps);
    }
}

// a = [1.0, 6.5, -2.3433, 7.34534, 1.2143, -3.11112, 2.1, 6.11, 3.3, 3.4, -4.1];

test "\t roots,\t real polynomial \t roots \n" {
    const eps = 1e-5;

    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);
        const n = 7;

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, n);
        const x: []C = try allocator.alloc(C, n - 1);

        a[0] = 1.0;
        a[1] = 1.2;
        a[2] = -3.2;
        a[3] = 2.1;
        a[4] = 3.3;
        a[5] = 3.4;
        a[6] = -4.1;

        try roots(R, allocator, a, x);

        try std.testing.expectEqual(x.len, 6);

        try std.testing.expectApproxEqAbs(@as(R, 0.5220415662921), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.3791318730909), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.6516102028100), x[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.4675374388097), x[3].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.5220415662921), x[4].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.6516102028100), x[5].re, eps);

        try std.testing.expectApproxEqAbs(@as(R, -0.4543195994722), x[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.7004994726793), x[2].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[3].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.4543195994722), x[4].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.7004994726793), x[5].im, eps);

        // print("\n", .{});
        // print("Complex Polynomial Residuals\n", .{});
        // var tmp = eval(R, C, a, x[0]); var hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
        // tmp = eval(R, C, a, x[1]); hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
        // tmp = eval(R, C, a, x[2]); hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
        // tmp = eval(R, C, a, x[3]); hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
        // tmp = eval(R, C, a, x[4]); hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
        // tmp = eval(R, C, a, x[5]); hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
    }
}

test "\t roots,\t real polynomial \t upperConvexHull \n" {
    const eps = 1e-5;

    inline for (.{ f32, f64 }) |R| {

        // const C: type = Complex(R);
        const n = 7;

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, n);
        var k: []usize = try allocator.alloc(usize, n);
        var hull: []R = try allocator.alloc(R, n);

        a[0] = 1.0;
        a[1] = 1.2;
        a[2] = -3.2;
        a[3] = 2.1;
        a[4] = 3.3;
        a[5] = 3.4;
        a[6] = -4.1;

        upperConvexHull(R, a, &k, &hull);

        try std.testing.expectEqual(@as(usize, 3), k.len);
        try std.testing.expectEqual(@as(usize, 0), k[0]);
        try std.testing.expectEqual(@as(usize, 2), k[1]);
        try std.testing.expectEqual(@as(usize, 6), k[2]);

        try std.testing.expectEqual(@as(usize, 3), hull.len);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), hull[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.678071905), hull[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, 2.035623909), hull[2], eps);
    }
}

test "\t roots,\t real polynomial \t computeModulii \n" {
    const eps = 1e-5;

    inline for (.{f64}) |R| {
        const n = 7;

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const a = try allocator.alloc(R, n);
        a[0] = 1.0;
        a[1] = 1.2;
        a[2] = -3.2;
        a[3] = 2.1;
        a[4] = 3.3;
        a[5] = 3.4;
        a[6] = -4.1;

        const k: []usize = try allocator.alloc(usize, 3);
        k[0] = 0;
        k[1] = 2;
        k[2] = 6;

        const u: []R = try allocator.alloc(R, 2);
        computeModulii(R, a, k, u);

        try std.testing.expectEqual(@as(usize, 2), u.len);
        try std.testing.expectApproxEqAbs(@as(R, 0.55901699437), u[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.939921384265), u[1], eps);
    }
}

test "\t roots,\t real polynomial \t initialRootEst \n" {
    const eps = 1e-5;

    inline for (.{f64}) |R| {
        const C: type = Complex(R);

        const n = 7;

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        const k: []usize = try allocator.alloc(usize, 3);
        k[0] = 0;
        k[1] = 2;
        k[2] = 6;

        const u: []R = try allocator.alloc(R, 2);
        u[0] = 0.55901699437;
        u[1] = 0.939921384265;

        const x: []C = try allocator.alloc(C, n - 1);
        initialRootEst(R, k, u, x);

        try std.testing.expectEqual(x.len, n - 1);
        try std.testing.expectApproxEqAbs(@as(R, 0.098100656), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.098100656), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.319821335), x[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.883836253), x[3].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.319821335), x[4].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.883836253), x[5].re, eps);

        try std.testing.expectApproxEqAbs(@as(R, -0.550341949), x[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.550341949), x[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.883836253), x[2].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.319821335), x[3].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.883836253), x[4].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.319821335), x[5].im, eps);
    }
}

test "\t roots,\t real polynomial \t aberthIterations \n" {
    const eps = 1e-5;

    inline for (.{f64}) |R| {
        const C: type = Complex(R);

        const n = 7;

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try allocator.alloc(R, n);
        a[0] = 1.0;
        a[1] = 1.2;
        a[2] = -3.2;
        a[3] = 2.1;
        a[4] = 3.3;
        a[5] = 3.4;
        a[6] = -4.1;

        const s_tilde = try allocator.alloc(R, n);
        for (a, 0..) |aval, i| {
            s_tilde[i] = fabs(R, aval) * @as(R, @floatFromInt((4 * i + 1)));
        }

        var stop_flags: []bool = try allocator.alloc(bool, n - 1);
        stop_flags[0] = false;
        stop_flags[1] = false;
        stop_flags[2] = false;
        stop_flags[3] = false;
        stop_flags[4] = false;
        stop_flags[5] = false;

        var x: []C = try allocator.alloc(C, n - 1);
        x[0].re = 0.098100656;
        x[1].re = -0.098100656;
        x[2].re = -0.319821335;
        x[3].re = 0.883836253;
        x[4].re = 0.319821335;
        x[5].re = -0.883836253;

        x[0].im = -0.550341949;
        x[1].im = 0.550341949;
        x[2].im = -0.883836253;
        x[3].im = -0.319821335;
        x[4].im = 0.883836253;
        x[5].im = 0.319821335;

        aberth(R, 11, a, s_tilde, stop_flags, x);

        try std.testing.expectEqual(x.len, n - 1);
        try std.testing.expectApproxEqAbs(@as(R, 0.52204156629), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.37913187309), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.65161020281), x[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.46753743881), x[3].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.52204156629), x[4].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.65161020281), x[5].re, eps);

        try std.testing.expectApproxEqAbs(@as(R, -0.45431959947), x[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.70049947268), x[2].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x[3].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.45431959947), x[4].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.70049947268), x[5].im, eps);
    }
}

// a = [1.0, 6.5, -2.3433, 7.34534, 1.2143, -3.11112, 2.1, 6.11, 3.3, 3.4, -4.1];

test "\t roots,\t complex polynomial \t roots \n" {
    const eps = 1e-5;

    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);
        const n = 7;

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try allocator.alloc(C, n);
        const x: []C = try allocator.alloc(C, n - 1);

        a[0].re = 1.0;
        a[1].re = 1.2;
        a[2].re = -3.2;
        a[3].re = 2.1;
        a[4].re = 3.3;
        a[5].re = 3.4;
        a[6].re = -4.1;

        a[0].im = -1.1;
        a[1].im = 1.12;
        a[2].im = 2.12;
        a[3].im = 1.0;
        a[4].im = 0.13;
        a[5].im = -1.1;
        a[6].im = -1.2;

        try roots(C, allocator, a, x);

        try std.testing.expectEqual(x.len, 6);

        try std.testing.expectApproxEqAbs(@as(R, 0.63240705162203), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.50361092865494), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.74808974891161), x[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.33794307164452), x[3].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.46969195349043), x[4].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.49683454987538), x[5].re, eps);

        try std.testing.expectApproxEqAbs(@as(R, -0.74783961161205), x[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.68097033792810), x[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.72824169200506), x[2].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.18140096660942), x[3].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.26056897926104), x[4].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.24525802153055), x[5].im, eps);

        // print("\n", .{});
        // print("Complex Polynomial Residuals\n", .{});
        // var tmp = eval(C, C, a, x[0]); var hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
        // tmp = eval(C, C, a, x[1]); hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
        // tmp = eval(C, C, a, x[2]); hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
        // tmp = eval(C, C, a, x[3]); hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
        // tmp = eval(C, C, a, x[4]); hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
        // tmp = eval(C, C, a, x[5]); hyp = fabs(C, tmp); print("{e:15.9}\n", .{hyp});
    }
}

test "\t roots,\t complex polynomial \t upperConvexHull \n" {
    const eps = 1e-5;

    inline for (.{ f32, f64 }) |R| {
        const C: type = Complex(R);
        const n = 7;

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try allocator.alloc(C, n);
        var k: []usize = try allocator.alloc(usize, n);
        var hull: []R = try allocator.alloc(R, n);

        a[0].re = 1.0;
        a[1].re = 1.2;
        a[2].re = -3.2;
        a[3].re = 2.1;
        a[4].re = 3.3;
        a[5].re = 3.4;
        a[6].re = -4.1;

        a[0].im = -1.1;
        a[1].im = 1.12;
        a[2].im = 2.12;
        a[3].im = 1.0;
        a[4].im = 0.13;
        a[5].im = -1.1;
        a[6].im = -1.2;

        upperConvexHull(C, a, &k, &hull);

        try std.testing.expectEqual(@as(usize, 3), k.len);
        try std.testing.expectEqual(@as(usize, 0), k[0]);
        try std.testing.expectEqual(@as(usize, 2), k[1]);
        try std.testing.expectEqual(@as(usize, 6), k[2]);

        try std.testing.expectEqual(@as(usize, 3), hull.len);
        try std.testing.expectApproxEqAbs(@as(R, 0.572023185), hull[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.940558204), hull[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, 2.094912279), hull[2], eps);
    }
}

test "\t roots,\t complex polynomial \t computeModulii \n" {
    const eps = 1e-5;

    inline for (.{f64}) |R| {
        const C: type = Complex(R);

        const n = 7;

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try allocator.alloc(C, n);
        a[0].re = 1.0;
        a[1].re = 1.2;
        a[2].re = -3.2;
        a[3].re = 2.1;
        a[4].re = 3.3;
        a[5].re = 3.4;
        a[6].re = -4.1;

        a[0].im = -1.1;
        a[1].im = 1.12;
        a[2].im = 2.12;
        a[3].im = 1.0;
        a[4].im = 0.13;
        a[5].im = -1.1;
        a[6].im = -1.2;

        var k: []usize = try allocator.alloc(usize, 3);
        k[0] = 0;
        k[1] = 2;
        k[2] = 6;

        const u: []R = try allocator.alloc(R, 2);
        computeModulii(C, a, k, u);

        try std.testing.expectEqual(@as(usize, 2), u.len);
        try std.testing.expectApproxEqAbs(@as(R, 0.62232171385), u[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.97360702387), u[1], eps);
    }
}

test "\t roots,\t complex polynomial \t initialRootEst \n" {
    const eps = 1e-5;

    inline for (.{f64}) |R| {
        const C: type = Complex(R);

        const n = 7;

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var k: []usize = try allocator.alloc(usize, 3);
        k[0] = 0;
        k[1] = 2;
        k[2] = 6;

        var u: []R = try allocator.alloc(R, 2);
        u[0] = 0.62232171385;
        u[1] = 0.97360702387;

        const x: []C = try allocator.alloc(C, n - 1);
        initialRootEst(R, k, u, x);

        try std.testing.expectEqual(@as(usize, n - 1), x.len);
        try std.testing.expectApproxEqAbs(@as(R, 0.109209861), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.109209861), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.331283342), x[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.915511869), x[3].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.331283342), x[4].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.915511869), x[5].re, eps);

        try std.testing.expectApproxEqAbs(@as(R, -0.612664281), x[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.612664281), x[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.915511869), x[2].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.331283342), x[3].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.915511869), x[4].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.331283342), x[5].im, eps);
    }
}

test "\t roots,\t complex polynomial \t aberthIterations \n" {
    const eps = 1e-5;

    inline for (.{f64}) |R| {
        const C: type = Complex(R);

        const n = 7;

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var a = try allocator.alloc(C, n);

        a[0].re = 1.0;
        a[1].re = 1.2;
        a[2].re = -3.2;
        a[3].re = 2.1;
        a[4].re = 3.3;
        a[5].re = 3.4;
        a[6].re = -4.1;

        a[0].im = -1.1;
        a[1].im = 1.12;
        a[2].im = 2.12;
        a[3].im = 1.0;
        a[4].im = 0.13;
        a[5].im = -1.1;
        a[6].im = -1.2;

        var s_tilde = try allocator.alloc(R, n);
        for (a, 0..) |aval, i| {
            s_tilde[i] = fabs(C, aval) * @as(R, @floatFromInt((4 * i + 1)));
        }

        var stop_flags: []bool = try allocator.alloc(bool, n - 1);
        stop_flags[0] = false;
        stop_flags[1] = false;
        stop_flags[2] = false;
        stop_flags[3] = false;
        stop_flags[4] = false;
        stop_flags[5] = false;

        var x: []C = try allocator.alloc(C, n - 1);
        x[0].re = 0.109209861;
        x[1].re = -0.109209861;
        x[2].re = -0.331283342;
        x[3].re = 0.915511869;
        x[4].re = 0.331283342;
        x[5].re = -0.915511869;

        x[0].im = -0.612664281;
        x[1].im = 0.612664281;
        x[2].im = -0.915511869;
        x[3].im = -0.331283342;
        x[4].im = 0.915511869;
        x[5].im = 0.331283342;

        aberth(C, 11, a, s_tilde, stop_flags, x);

        try std.testing.expectApproxEqAbs(@as(R, 0.63240705162203), x[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.50361092865494), x[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.74808974891162), x[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.33794307164455), x[3].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.46969195349043), x[4].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.49683454987538), x[5].re, eps);

        try std.testing.expectApproxEqAbs(@as(R, -0.74783961161205), x[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.68097033792810), x[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.72824169200506), x[2].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.18140096660942), x[3].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.26056897926104), x[4].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.24525802153055), x[5].im, eps);
    }
}
