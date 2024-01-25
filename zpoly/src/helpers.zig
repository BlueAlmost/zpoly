const std = @import("std");
const print = std.debug.print;
const math = std.math;
const Complex = std.math.complex.Complex;

const zpoly = @import("../zpoly.zig");

// Helper comptime type functions -----------------------------------------------------

pub fn EitherCmpx(comptime A: type, comptime B: type) type {
    if (A == B) {
        return A;
    } else if (A == Complex(B)) {
        return A;
    } else if (B == Complex(A)) {
        return B;
    } else {
        @compileError("incompatible real/complex pairings");
    }
}

pub fn ToCmpx(comptime T: type) type {
    switch (T) {
        f32, f64 => {
            return Complex(T);
        },
        Complex(f32), Complex(f64) => {
            return T;
        },
        else => {
            @compileError("unexpected type");
        },
    }
}

pub fn ValueType(comptime T: type) type {
    switch (T) {
        f32, f64 => {
            return T;
        },
        Complex(f32) => {
            return f32;
        },
        Complex(f64) => {
            return f64;
        },
        else => {
            @compileError("unexpected type");
        },
    }
}
