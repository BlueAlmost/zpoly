const std = @import("std");
const math = std.math;

// helper functions
pub const EitherCmpx = @import("src/helpers.zig").EitherCmpx;
pub const ToCmpx = @import("src/helpers.zig").ToCmpx;
pub const ValueType = @import("src/helpers.zig").ValueType;

// complex mathematics
pub const cmath = @import("src/complex_math.zig");

// complex number formatting
pub const fmtCi = @import("src/complex_fmt.zig").fmtCi;
pub const fmtCp = @import("src/complex_fmt.zig").fmtCp;

// polynomial functions
pub const prt_descend = @import("src/polynomial.zig").prt_descend;
pub const prt_ascend = @import("src/polynomial.zig").prt_ascend;
pub const monic = @import("src/polynomial.zig").monic;

// polynomial mathematical functions
pub const fabs = @import("src/math.zig").fabs;
pub const add = @import("src/math.zig").add;
pub const mult = @import("src/math.zig").conv;

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
pub const wrapper_roots = @import("src/roots.zig").wrapper_roots;
