const std = @import("std");
const math = std.math;

// container
pub const Polynomial = @import("src/polynomial.zig").Polynomial;

// helper functions
pub const EitherCmpx = @import("src/helpers.zig").EitherCmpx;
pub const ToCmpx = @import("src/helpers.zig").ToCmpx;
pub const ValueType = @import("src/helpers.zig").ValueType;

// polynomial mathematical functions
pub const add = @import("src/math.zig").add;
pub const mult = @import("src/math.zig").mult;

pub const eval = @import("src/math.zig").eval;
pub const evalRev = @import("src/math.zig").evalRev;
pub const evalDeriv = @import("src/math.zig").evalDeriv;
pub const evalDerivRev = @import("src/math.zig").evalDerivRev;

// root functions
pub const upperConvexHull = @import("src/roots.zig").upperConvexHull;
pub const computeModulii = @import("src/roots.zig").computeModulii;
pub const initialRootEst = @import("src/roots.zig").initialRootEst;
pub const aberth = @import("src/roots.zig").aberth;
pub const roots = @import("src/roots.zig").roots;
