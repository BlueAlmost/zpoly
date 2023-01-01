const std = @import("std");
const print = std.debug.print;
const math = std.math;
const Complex = std.math.complex.Complex;

const Allocator = std.mem.Allocator;

// NOTE:  many loops below use the "order" field, because we only care
// about coefficients upto the leading coefficient of the polynomial
// (i.e. a[0] up to a[p] where "p" is the order of our polynomial)

pub fn Polynomial(comptime T: type) type {
    comptime var U: type = usize;
    comptime var R: type = undefined;

    switch (T) {
        f32 => {
            R = f32;
        },
        f64 => {
            R = f64;
        },
        Complex(f32) => {
            R = f32;
        },
        Complex(f64) => {
            R = f64;
        },
        else => @compileError("requested Polynomial type is not allowed"),
    }

    return struct {
        const Self = @This();
        maxlen: U,
        val: []T,

        pub fn init(allocator: Allocator, maxlen: U, len: U) !Self {
            if (maxlen < len) {
                return error.Order_Exceeds_Maxlen;
            }

            var val = try allocator.alloc(T, maxlen);
            val.len = len;

            return Self{
                .maxlen = maxlen,
                .val = val,
            };
        }

        pub fn fill(self: Self, new_val: T) void {
            var i: usize = 0;
            while (i < self.val.len) : (i += 1) {
                self.val[i] = new_val;
            }
        }

        pub fn monic(self: *Self, eps: R) void {
            switch (T) {
                f32, f64 => {
                    var i: usize = self.val.len - 1;
                    while (i > 0) : (i -= 1) {
                        if (@fabs(self.val[i]) < eps) {
                            self.val.len -= 1;
                        }
                    }

                    var alpha: T = 1 / self.val[self.val.len - 1];
                    i = 0;
                    while (i < self.val.len - 1) : (i += 1) {
                        self.val[i] *= alpha;
                    }
                    self.val[self.val.len - 1] = 1;
                },

                Complex(f32), Complex(f64) => {
                    var i: usize = self.val.len - 1;
                    while (i > 0) : (i -= 1) {
                        if (math.complex.abs(self.val[i]) < eps) {
                            self.val.len -= 1;
                        }
                    }

                    var val: T = self.val[self.val.len - 1];
                    var beta: R = val.re * val.re + val.im * val.im;
                    var alpha = T.init(val.re / beta, -val.im / beta);

                    i = 0;
                    while (i < self.val.len - 1) : (i += 1) {
                        var tmp = self.val[i];
                        self.val[i].re = alpha.re * tmp.re - alpha.im * tmp.im;
                        self.val[i].im = alpha.im * tmp.re + alpha.re * tmp.im;
                    }
                    self.val[self.val.len - 1].re = 1;
                    self.val[self.val.len - 1].im = 0;
                },
                else => {
                    @compileError("unknown Polynomial type");
                },
            }
        }

        pub fn neg(self: Self) void {
            switch (T) {
                f32, f64 => {
                    var i: usize = 0;
                    while (i < self.val.len) : (i += 1) {
                        self.val[i] = -self.val[i];
                    }
                },

                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    while (i < self.val.len) : (i += 1) {
                        self.val[i].re = -self.val[i].re;
                        self.val[i].im = -self.val[i].im;
                    }
                },
                else => @compileError("requested Polynomial type is not allowed"),
            }
        }

        pub fn ones(self: Self) void {
            switch (T) {
                f32, f64 => {
                    var i: usize = 0;
                    while (i < self.val.len) : (i += 1) {
                        self.val[i] = 1;
                    }
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    while (i < self.val.len) : (i += 1) {
                        self.val[i].re = 1;
                        self.val[i].im = 0;
                    }
                },
                else => @compileError("requested Polynomial type is not allowed"),
            }
        }

        pub fn prt(self: *Self) void {
            switch (T) {
                f32, f64 => {
                    // print representation as
                    //   a_0 + a_1*x + a_2*x]^2 + ... a_n*x^n

                    print(" \n", .{});
                    var i: usize = 0;
                    while (i < self.val.len) : (i += 1) {
                        var val: T = self.val[i];

                        if (val < 0) {
                            print(" {d:.4} x^{d}", .{ val, i });
                        } else if ((val >= 0) and (i == 0)) {
                            print("   {d:.4} x^{d}", .{ val, i });
                        } else if (val >= 0) {
                            print(" + {d:.4} x^{d}", .{ val, i });
                        }
                    }
                    print(" \n", .{});
                },

                Complex(f32), Complex(f64) => {
                    print(" \n", .{});

                    var i: usize = 0;
                    while (i < self.val.len) : (i += 1) {
                        var val: T = self.val[i];

                        if (i == 0) {
                            print(" ({d:.4}, {d:.4}) x^{d}", .{ val.re, val.im, i });
                        } else {
                            print(" + ({d:.4}, {d:.4}) x^{d}", .{ val.re, val.im, i });
                        }
                    }
                    print(" \n", .{});
                },
                else => @compileError("requested Polynomial type is not allowed"),
            }
        }

        pub fn scale(self: Self, alpha: anytype) void {
            comptime {
                if ((@TypeOf(alpha) != R) and (@TypeOf(alpha) != comptime_float)) {
                    @compileError("scaling parameter must be real");
                }
            }

            switch (T) {
                f32, f64 => {
                    var i: usize = 0;
                    while (i < self.val.len) : (i += 1) {
                        self.val[i] *= alpha;
                    }
                },

                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    while (i < self.val.len) : (i += 1) {
                        self.val[i].re *= alpha;
                        self.val[i].im *= alpha;
                    }
                },
                else => @compileError("requested Polynomial type is not allowed"),
            }
        }

        pub fn set_len(self: Self, len: usize) void {
            if (self.maxlen < len) {
                return error.Maxlen_Cannot_Support_Requested_Polynomial_Length;
            } else {
                self.val.len = len;
            }
        }

        pub fn trim(self: *Self, eps: R) void {
            switch (T) {
                f32, f64 => {
                    var i: usize = self.val.len - 1;

                    while (i > 0) : (i -= 1) {
                        if (@fabs(self.val[i]) < eps) {
                            self.val.len -= 1;
                        }
                    }
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = self.val.len - 1;

                    while (i > 0) : (i -= 1) {
                        if (math.complex.abs(self.val[i]) < eps) {
                            self.val.len -= 1;
                        }
                    }
                },
                else => {
                    @compileError("unknown Polynomial type");
                },
            }
        }

        pub fn zeros(self: Self) void {
            switch (T) {
                f32, f64 => {
                    var i: usize = 0;
                    while (i < self.val.len) : (i += 1) {
                        self.val[i] = 0;
                    }
                },
                Complex(f32), Complex(f64) => {
                    var i: usize = 0;
                    while (i < self.val.len) : (i += 1) {
                        self.val[i].re = 0;
                        self.val[i].im = 0;
                    }
                },
                else => @compileError("requested Polynomial type is not allowed"),
            }
        }
    };
}

//------------- TEST ------------------------------------------

test "\t polynomial methods,\t fill - real \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(R).init(allocator, 5, 5);
        x.fill(12.34);

        try std.testing.expectApproxEqAbs(@as(R, 12.34), x.val[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 12.34), x.val[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, 12.34), x.val[2], eps);
        try std.testing.expectApproxEqAbs(@as(R, 12.34), x.val[3], eps);
        try std.testing.expectApproxEqAbs(@as(R, 12.34), x.val[4], eps);
    }
}

test "\t polynomial methods,\t fill - complex \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(C).init(allocator, 3, 3);
        x.fill(C.init(1.23, 4.56));

        try std.testing.expectApproxEqAbs(@as(R, 1.23), x.val[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.56), x.val[0].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.23), x.val[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.56), x.val[1].im, eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.23), x.val[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.56), x.val[2].im, eps);
    }
}

test "\t polynomial methods,\t monic - real \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(R).init(allocator, 6, 5);
        x.val[0] = 2.1;
        x.val[1] = 4.2;
        x.val[2] = -0.21;
        x.val[3] = 2.1;
        x.val[4] = 0.0001;

        x.monic(0.001);

        try std.testing.expectEqual(x.val.len, 4);
        try std.testing.expectApproxEqAbs(@as(R, 1.0), x.val[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 2.0), x.val[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, -0.1), x.val[2], eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.0), x.val[3], eps);
    }
}

test "\t polynomial methods,\t monic - complex \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(C).init(allocator, 6, 5);
        x.val[0].re = 2.1;
        x.val[0].im = 2.1;

        x.val[1].re = 4.2;
        x.val[1].im = -21.0;

        x.val[2].re = 0.21;
        x.val[2].im = 0.21;

        x.val[3].re = 2.1;
        x.val[3].im = -3.3;

        x.val[4].re = 0.0001;
        x.val[4].im = 0.0001;

        x.monic(0.001);

        try std.testing.expectEqual(@as(usize, 4), x.val.len);
        try std.testing.expectApproxEqAbs(@as(R, -0.164705882353), x.val[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.741176405888), x.val[0].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 5.105882352941), x.val[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -1.976470588235), x.val[1].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, -0.016470588235), x.val[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.074117647058), x.val[2].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 1.0), x.val[3].re, eps);
    }
}

test "\t polynomial methods,\t neg - real \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(R).init(allocator, 2, 2);
        x.val[0] = 1.23;
        x.val[1] = 4.56;
        x.neg();

        try std.testing.expectApproxEqAbs(@as(R, -1.23), x.val[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, -4.56), x.val[1], eps);
    }
}

test "\t polynomial methods,\t neg - complex \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(C).init(allocator, 2, 2);
        x.val[0] = C.init(1.23, 4.56);
        x.val[1] = C.init(7.89, -4.56);
        x.neg();

        try std.testing.expectApproxEqAbs(@as(R, -1.23), x.val[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -4.56), x.val[0].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, -7.89), x.val[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 4.56), x.val[1].im, eps);
    }
}

test "\t polynomial methods,\t ones - real \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(R).init(allocator, 2, 2);
        x.ones();

        try std.testing.expectApproxEqAbs(@as(R, 1.0), x.val[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 1.0), x.val[1], eps);
    }
}

test "\t polynomial methods,\t ones - complex \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        const C = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(C).init(allocator, 2, 2);
        x.ones();

        try std.testing.expectApproxEqAbs(@as(R, 1.0), x.val[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x.val[0].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 1.0), x.val[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x.val[1].im, eps);
    }
}

test "\t polynomial methods,\t scale - real \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(R).init(allocator, 2, 2);
        x.val[0] = 1.1;
        x.val[1] = 3.3;

        var alpha: R = 2.0;

        x.scale(alpha);

        try std.testing.expectApproxEqAbs(@as(R, 2.2), x.val[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 6.6), x.val[1], eps);
    }
}

test "\t polynomial methods,\t scale - complex \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        comptime var C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(C).init(allocator, 2, 2);
        x.val[0] = C.init(1.1, 3.3);
        x.val[1] = C.init(1.1, -3.3);

        // var alpha: R = 2.0;
        // x.scale(alpha);
        x.scale(2.0);

        try std.testing.expectApproxEqAbs(@as(R, 2.2), x.val[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 6.6), x.val[0].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 2.2), x.val[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, -6.6), x.val[1].im, eps);
    }
}

test "\t polynomial methods,\t trim - real \n" {
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(R).init(allocator, 5, 5);
        x.fill(12.34);
        x.val[4] = 0.0001;
        x.val[3] = 0.0001;

        x.trim(0.01);

        try std.testing.expectEqual(@as(usize, 3), x.val.len);
    }
}

test "\t polynomial methods,\t trim - complex \n" {
    inline for (.{ Complex(f32), Complex(f64) }) |C| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(C).init(allocator, 5, 5);
        x.fill(C.init(12.34, 56.78));
        x.val[4].re = 0.000001;
        x.val[4].im = 0.000001;
        x.val[3].re = 0.000001;
        x.val[3].im = 0.000001;

        x.trim(0.01);

        try std.testing.expectEqual(@as(usize, 3), x.val.len);
    }
}

test "\t polynomial methods,\t zeros - real \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(R).init(allocator, 3, 3);
        x.zeros();

        try std.testing.expectApproxEqAbs(@as(R, 0.0), x.val[0], eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x.val[1], eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x.val[2], eps);
    }
}

test "\t polynomial methods,\t zeros - complex \n" {
    const eps = 1e-6;
    inline for (.{ f32, f64 }) |R| {
        comptime var C: type = Complex(R);

        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();

        var x = try Polynomial(C).init(allocator, 3, 3);
        x.zeros();

        try std.testing.expectApproxEqAbs(@as(R, 0.0), x.val[0].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x.val[0].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 0.0), x.val[1].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x.val[1].im, eps);

        try std.testing.expectApproxEqAbs(@as(R, 0.0), x.val[2].re, eps);
        try std.testing.expectApproxEqAbs(@as(R, 0.0), x.val[2].im, eps);
    }
}
