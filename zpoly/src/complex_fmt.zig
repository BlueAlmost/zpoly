const std = @import("std");
const Formatter = std.fmt.Formatter;
const formatFloat = std.fmt.formatFloat;
const FormatOptions = std.fmt.FormatOptions;
const invalidFmtError = std.fmt.invalidFmtError;

const Complex = std.math.complex.Complex;

pub fn fmtCi(data: anytype) Formatter(formatCi) {
    return .{.data = Complex(f64).init(@as(f64, @floatCast(data.re)),
        @as(f64, @floatCast(data.im)))};
}

pub fn formatCi( data: Complex(f64), comptime fmt: []const u8, options: FormatOptions, writer: anytype) !void {

    // format complex number using the "i" notation

    var buf: [128]u8 = undefined;

    const prefix_re: u8 = if (data.re > 0) ' ' else '-';
    const prefix_im: u8 = if (data.im > 0) '+' else '-';

    const abs_data_re: @TypeOf(data.re) = @abs(data.re);
    const abs_data_im: @TypeOf(data.im) = @abs(data.im);

    _ = try writer.writeByte(prefix_re);
    if (fmt.len == 0 or comptime std.mem.eql(u8, fmt, "e")) {

        var slice = try formatFloat(&buf, abs_data_re, .{.mode = .scientific, .precision = options.precision});
        _ = try writer.write(slice);

        _ = try writer.writeByte(prefix_im);
        slice = try formatFloat(&buf, abs_data_im, .{.mode = .scientific, .precision = options.precision});
        _ = try writer.write(slice);
    
    } else if (comptime std.mem.eql(u8, fmt, "d")) {

        var slice = try formatFloat(&buf, abs_data_re, .{.mode = .decimal, .precision = options.precision});
        _ = try writer.write(slice);

        _ = try writer.writeByte(prefix_im);
        slice = try formatFloat(&buf, abs_data_im, .{.mode = .decimal, .precision = options.precision});
        _ = try writer.write(slice);
    } else {
        invalidFmtError(fmt, data);
    }

    try writer.writeAll("i");
}


pub fn fmtCp(data: anytype) Formatter(formatCp) {
    return .{.data = Complex(f64).init(@as(f64, @floatCast(data.re)),
        @as(f64, @floatCast(data.im)))};
}

pub fn formatCp( data: Complex(f64), comptime fmt: []const u8, options: FormatOptions, writer: anytype) !void {

    // format complex number using the "i" notation

    var buf: [128]u8 = undefined;

    const prefix_re: u8 = if (data.re > 0) ' ' else '-';
    const prefix_im: u8 = if (data.im > 0) ' ' else '-';

    const abs_data_re: @TypeOf(data.re) = @abs(data.re);
    const abs_data_im: @TypeOf(data.im) = @abs(data.im);

    try writer.writeAll("(");

    _ = try writer.writeByte(prefix_re);
    if (fmt.len == 0 or comptime std.mem.eql(u8, fmt, "e")) {

        var slice = try formatFloat(&buf, abs_data_re, .{.mode = .scientific, .precision = options.precision});
        _ = try writer.write(slice);
    
        try writer.writeAll(",");

        _ = try writer.writeByte(prefix_im);
        slice = try formatFloat(&buf, abs_data_im, .{.mode = .scientific, .precision = options.precision});
        _ = try writer.write(slice);
    
    } else if (comptime std.mem.eql(u8, fmt, "d")) {

        var slice = try formatFloat(&buf, abs_data_re, .{.mode = .decimal, .precision = options.precision});
        _ = try writer.write(slice);

        try writer.writeAll(",");

        _ = try writer.writeByte(prefix_im);
        slice = try formatFloat(&buf, abs_data_im, .{.mode = .decimal, .precision = options.precision});
        _ = try writer.write(slice);
    } else {
        invalidFmtError(fmt, data);
    }

    try writer.writeAll(")");
}


test "complex i decimal format test\n" {
    var buf: [100]u8 = undefined;

    const nn32 = Complex(f32).init(-11.111, -222.222);
    const np32 = Complex(f32).init(-11.111, 222.222);
    const pn32 = Complex(f32).init(11.111, -222.222);
    const pp32 = Complex(f32).init(11.111, 222.222);

    var slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(nn32)});
    try std.testing.expect(std.mem.eql(u8, "-11.111-222.222i", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(np32)});
    try std.testing.expect(std.mem.eql(u8, "-11.111+222.222i", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(pn32)});
    try std.testing.expect(std.mem.eql(u8, " 11.111-222.222i", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(pp32)});
    try std.testing.expect(std.mem.eql(u8, " 11.111+222.222i", slice));

    const nn64 = Complex(f64).init(-11.111, -222.222);
    const np64 = Complex(f64).init(-11.111, 222.222);
    const pn64 = Complex(f64).init(11.111, -222.222);
    const pp64 = Complex(f64).init(11.111, 222.222);

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(nn64)});
    try std.testing.expect(std.mem.eql(u8, "-11.111-222.222i", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(np64)});
    try std.testing.expect(std.mem.eql(u8, "-11.111+222.222i", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(pn64)});
    try std.testing.expect(std.mem.eql(u8, " 11.111-222.222i", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCi(pp64)});
    try std.testing.expect(std.mem.eql(u8, " 11.111+222.222i", slice));
}

test "complex p exponential format test\n" {
    var buf: [100]u8 = undefined;

    const nn32 = Complex(f32).init(-11.111, -222.222);
    const np32 = Complex(f32).init(-11.111, 222.222);
    const pn32 = Complex(f32).init(11.111, -222.222);
    const pp32 = Complex(f32).init(11.111, 222.222);

    var slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(nn32)});
    try std.testing.expect(std.mem.eql(u8, "(-1.111e1,-2.222e2)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(np32)});
    try std.testing.expect(std.mem.eql(u8, "(-1.111e1, 2.222e2)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(pn32)});
    try std.testing.expect(std.mem.eql(u8, "( 1.111e1,-2.222e2)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(pp32)});
    try std.testing.expect(std.mem.eql(u8, "( 1.111e1, 2.222e2)", slice));

    const nn64 = Complex(f64).init(-11.111, -222.222);
    const np64 = Complex(f64).init(-11.111, 222.222);
    const pn64 = Complex(f64).init(11.111, -222.222);
    const pp64 = Complex(f64).init(11.111, 222.222);

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(nn64)});
    try std.testing.expect(std.mem.eql(u8, "(-1.111e1,-2.222e2)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(np64)});
    try std.testing.expect(std.mem.eql(u8, "(-1.111e1, 2.222e2)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(pn64)});
    try std.testing.expect(std.mem.eql(u8, "( 1.111e1,-2.222e2)", slice));

    slice = try std.fmt.bufPrint(&buf, "{e:6.3}", .{fmtCp(pp64)});
    try std.testing.expect(std.mem.eql(u8, "( 1.111e1, 2.222e2)", slice));
}

test "complex p decimal format test\n" {
    var buf: [100]u8 = undefined;

    const nn32 = Complex(f32).init(-11.111, -222.222);
    const np32 = Complex(f32).init(-11.111, 222.222);
    const pn32 = Complex(f32).init(11.111, -222.222);
    const pp32 = Complex(f32).init(11.111, 222.222);

    var slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(nn32)});
    try std.testing.expect(std.mem.eql(u8, "(-11.111,-222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(np32)});
    try std.testing.expect(std.mem.eql(u8, "(-11.111, 222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(pn32)});
    try std.testing.expect(std.mem.eql(u8, "( 11.111,-222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(pp32)});
    try std.testing.expect(std.mem.eql(u8, "( 11.111, 222.222)", slice));

    const nn64 = Complex(f64).init(-11.111, -222.222);
    const np64 = Complex(f64).init(-11.111, 222.222);
    const pn64 = Complex(f64).init(11.111, -222.222);
    const pp64 = Complex(f64).init(11.111, 222.222);

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(nn64)});
    try std.testing.expect(std.mem.eql(u8, "(-11.111,-222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(np64)});
    try std.testing.expect(std.mem.eql(u8, "(-11.111, 222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(pn64)});
    try std.testing.expect(std.mem.eql(u8, "( 11.111,-222.222)", slice));

    slice = try std.fmt.bufPrint(&buf, "{d:.3}", .{fmtCp(pp64)});
    try std.testing.expect(std.mem.eql(u8, "( 11.111, 222.222)", slice));
}
