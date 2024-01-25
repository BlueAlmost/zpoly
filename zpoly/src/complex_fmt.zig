//
// custom formatters for printing complex numbers
//

const std = @import("std");

const Formatter = std.fmt.Formatter;
const Complex = std.math.complex.Complex;

const formatFloatScientific = std.fmt.formatFloatScientific;
const formatFloatDecimal = std.fmt.formatFloatDecimal;
const formatFloatHexadecimal = std.fmt.formatFloatHexadecimal;
const invalidFmtErr = std.fmt.invalidFmtErr;
const formatBuf = std.fmt.formatBuf;
const print = std.debug.print;

pub fn fmtCi(data: anytype) Formatter(formatCi) {
    // format complex number using "i" notation (without parentheses)
    return .{ .data = Complex(f64).init(@as(f64, @floatCast(data.re)), @as(f64, @floatCast(data.im))) };
}

fn formatCi(data: Complex(f64), comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {

    // this buffer should be enough to display all decimal places of two decimal f64 numbers.
    var buf: [128]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    var buf_writer = fbs.writer();

    const prefix_re: u8 = if (data.re > 0) ' ' else '-';
    const prefix_im: u8 = if (data.im > 0) '+' else '-';

    const abs_data_re = @abs(data.re);
    const abs_data_im = @abs(data.im);

    try buf_writer.print("{c}", .{prefix_re});
    if (fmt.len == 0 or comptime std.mem.eql(u8, fmt, "e")) {
        try formatFloatScientific(abs_data_re, options, buf_writer);
        try buf_writer.print("{c}i*", .{prefix_im});
        try formatFloatScientific(abs_data_im, options, buf_writer);
    } else if (comptime std.mem.eql(u8, fmt, "d")) {
        try formatFloatDecimal(abs_data_re, options, buf_writer);
        try buf_writer.print("{c}i*", .{prefix_im});
        try formatFloatDecimal(abs_data_im, options, buf_writer);
    } else {
        invalidFmtErr(fmt, data);
    }
    return formatBuf(fbs.getWritten(), options, writer);
}

pub fn fmtCp(data: anytype) Formatter(formatCp) {
    // format complex number using parentheses notation ("p")
    return .{ .data = Complex(f64).init(@as(f64, @floatCast(data.re)), @as(f64, @floatCast(data.im))) };
}

fn formatCp(data: Complex(f64), comptime fmt: []const u8, options: std.fmt.FormatOptions, writer: anytype) !void {

    // this buffer should be enough to display all decimal places of two decimal f64 numbers.
    var buf: [128]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    var buf_writer = fbs.writer();

    const prefix_re: u8 = if (data.re > 0) ' ' else '-';
    const prefix_im: u8 = if (data.im > 0) ' ' else '-';

    const abs_data_re = @abs(data.re);
    const abs_data_im = @abs(data.im);

    try buf_writer.print("({c}", .{prefix_re});
    if (fmt.len == 0 or comptime std.mem.eql(u8, fmt, "e")) {
        try formatFloatScientific(abs_data_re, options, buf_writer);
        try buf_writer.print(", {c}", .{prefix_im});
        try formatFloatScientific(abs_data_im, options, buf_writer);
    } else if (comptime std.mem.eql(u8, fmt, "d")) {
        try formatFloatDecimal(abs_data_re, options, buf_writer);
        try buf_writer.print(", {c}", .{prefix_im});
        try formatFloatDecimal(abs_data_im, options, buf_writer);
    } else {
        invalidFmtErr(fmt, data);
    }
    try buf_writer.print(")", .{});
    return formatBuf(fbs.getWritten(), options, writer);
}

test "complex i exponential format test\n" {
    var buf: [100]u8 = undefined;

    const nn32 = Complex(f32).init(-11.111, -222.222);
    const np32 = Complex(f32).init(-11.111, 222.222);
    const pn32 = Complex(f32).init(11.111, -222.222);
    const pp32 = Complex(f32).init(11.111, 222.222);

    var slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCi(nn32)});
    try std.testing.expect(std.mem.eql(u8, "-1.111e+01-i*2.222e+02", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCi(np32)});
    try std.testing.expect(std.mem.eql(u8, "-1.111e+01+i*2.222e+02", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCi(pn32)});
    try std.testing.expect(std.mem.eql(u8, " 1.111e+01-i*2.222e+02", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCi(pp32)});
    try std.testing.expect(std.mem.eql(u8, " 1.111e+01+i*2.222e+02", slice));

    const nn64 = Complex(f64).init(-11.111, -222.222);
    const np64 = Complex(f64).init(-11.111, 222.222);
    const pn64 = Complex(f64).init(11.111, -222.222);
    const pp64 = Complex(f64).init(11.111, 222.222);

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCi(nn64)});
    try std.testing.expect(std.mem.eql(u8, "-1.111e+01-i*2.222e+02", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCi(np64)});
    try std.testing.expect(std.mem.eql(u8, "-1.111e+01+i*2.222e+02", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCi(pn64)});
    try std.testing.expect(std.mem.eql(u8, " 1.111e+01-i*2.222e+02", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCi(pp64)});
    try std.testing.expect(std.mem.eql(u8, " 1.111e+01+i*2.222e+02", slice));
}

test "complex i decimal format test\n" {
    var buf: [100]u8 = undefined;

    const nn32 = Complex(f32).init(-11.111, -222.222);
    const np32 = Complex(f32).init(-11.111, 222.222);
    const pn32 = Complex(f32).init(11.111, -222.222);
    const pp32 = Complex(f32).init(11.111, 222.222);

    var slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(nn32)});
    try std.testing.expect(std.mem.eql(u8, "-11.111-i*222.222", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(np32)});
    try std.testing.expect(std.mem.eql(u8, "-11.111+i*222.222", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(pn32)});
    try std.testing.expect(std.mem.eql(u8, " 11.111-i*222.222", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(pp32)});
    try std.testing.expect(std.mem.eql(u8, " 11.111+i*222.222", slice));

    const nn64 = Complex(f64).init(-11.111, -222.222);
    const np64 = Complex(f64).init(-11.111, 222.222);
    const pn64 = Complex(f64).init(11.111, -222.222);
    const pp64 = Complex(f64).init(11.111, 222.222);

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(nn64)});
    try std.testing.expect(std.mem.eql(u8, "-11.111-i*222.222", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(np64)});
    try std.testing.expect(std.mem.eql(u8, "-11.111+i*222.222", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(pn64)});
    try std.testing.expect(std.mem.eql(u8, " 11.111-i*222.222", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(pp64)});
    try std.testing.expect(std.mem.eql(u8, " 11.111+i*222.222", slice));
}

test "complex p exponential format test\n" {
    var buf: [100]u8 = undefined;

    const nn32 = Complex(f32).init(-11.111, -222.222);
    const np32 = Complex(f32).init(-11.111, 222.222);
    const pn32 = Complex(f32).init(11.111, -222.222);
    const pp32 = Complex(f32).init(11.111, 222.222);

    var slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(nn32)});
    try std.testing.expect(std.mem.eql(u8, "(-1.111e+01, -2.222e+02)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(np32)});
    try std.testing.expect(std.mem.eql(u8, "(-1.111e+01,  2.222e+02)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(pn32)});
    try std.testing.expect(std.mem.eql(u8, "( 1.111e+01, -2.222e+02)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(pp32)});
    try std.testing.expect(std.mem.eql(u8, "( 1.111e+01,  2.222e+02)", slice));

    const nn64 = Complex(f64).init(-11.111, -222.222);
    const np64 = Complex(f64).init(-11.111, 222.222);
    const pn64 = Complex(f64).init(11.111, -222.222);
    const pp64 = Complex(f64).init(11.111, 222.222);

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(nn64)});
    try std.testing.expect(std.mem.eql(u8, "(-1.111e+01, -2.222e+02)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(np64)});
    try std.testing.expect(std.mem.eql(u8, "(-1.111e+01,  2.222e+02)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(pn64)});
    try std.testing.expect(std.mem.eql(u8, "( 1.111e+01, -2.222e+02)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(pp64)});
    try std.testing.expect(std.mem.eql(u8, "( 1.111e+01,  2.222e+02)", slice));
}

test "complex p decimal format test\n" {
    var buf: [100]u8 = undefined;

    const nn32 = Complex(f32).init(-11.111, -222.222);
    const np32 = Complex(f32).init(-11.111, 222.222);
    const pn32 = Complex(f32).init(11.111, -222.222);
    const pp32 = Complex(f32).init(11.111, 222.222);

    var slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(nn32)});
    try std.testing.expect(std.mem.eql(u8, "(-11.111, -222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(np32)});
    try std.testing.expect(std.mem.eql(u8, "(-11.111,  222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(pn32)});
    try std.testing.expect(std.mem.eql(u8, "( 11.111, -222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(pp32)});
    try std.testing.expect(std.mem.eql(u8, "( 11.111,  222.222)", slice));

    const nn64 = Complex(f64).init(-11.111, -222.222);
    const np64 = Complex(f64).init(-11.111, 222.222);
    const pn64 = Complex(f64).init(11.111, -222.222);
    const pp64 = Complex(f64).init(11.111, 222.222);

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(nn64)});
    try std.testing.expect(std.mem.eql(u8, "(-11.111, -222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(np64)});
    try std.testing.expect(std.mem.eql(u8, "(-11.111,  222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(pn64)});
    try std.testing.expect(std.mem.eql(u8, "( 11.111, -222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(pp64)});
    try std.testing.expect(std.mem.eql(u8, "( 11.111,  222.222)", slice));
}
