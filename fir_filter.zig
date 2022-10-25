//!FIR滤波器

const std = @import("std");

const WindowFunction = @import("window_function.zig");

///滤波
///<param name="b">滤波器分子系数向量</param>
///<param name="a">滤波器分母系数向量</param>
///<param name="x">N维输入数组</param>
///<returns></returns>
pub fn Filter(b: []f64, a: []f64, x: []f64) ![]f64 {
    const allocator = std.heap.page_allocator;

    var y = try allocator.alloc(f64, x.len);
    if (a[0] == 0) return y;
    var i: usize = undefined;
    if (a[0] != 1) {
        i = 1;
        while (i < b.len) : (i += 1) {
            b[i] /= a[0];
        }
        i = 1;
        while (i < a.len) : (i += 1) {
            a[i] /= a[0];
        }
        a[0] = 1;
    }
    y[0] = x[0] * b[0];
    i = 1;
    while (i < x.len) : (i += 1) {
        y[i] = 0;
        var j: usize = 0;
        while (j < b.len) : (j += 1) {
            if (i < j) break;
            y[i] += b[j] * x[i - j];
        }
        j = 1;
        while (j < a.len) : (j += 1) {
            if (i < j) break;
            y[i] -= a[j] * x[i - j];
        }
    }
    return y;
}

pub fn Lowpass(numtaps: usize, cutoff: f64, width: f64, fs: f64) ![]f64 {
    const allocator = std.heap.page_allocator;

    const atten = kaiser_atten(numtaps, width / (0.5 * fs));
    const beta = kaiser_beta(atten);
    const cutoff_ = cutoff / (0.5 * fs);
    const wc = std.math.pi * cutoff_;

    const hd = try ideal_lp(numtaps, wc);

    var hn = try allocator.alloc(f64, numtaps);
    const wn = try WindowFunction.kaiser(allocator, numtaps, beta);
    var i: usize = 0;
    while (i < numtaps) : (i += 1) {
        hn[i] = hd[i] * wn[i];
    }
    return hn;
}

pub fn Highpass(numtaps: usize, cutoff: f64, width: f64, fs: f64) ![]f64 {
    const allocator = std.heap.page_allocator;

    const atten = kaiser_atten(numtaps, width / (0.5 * fs));
    const beta = kaiser_beta(atten);
    const cutoff_ = cutoff / (0.5 * fs);
    const wc = std.math.pi * cutoff_;

    const hd1 = try ideal_lp(numtaps, std.math.pi);
    const hd2 = try ideal_lp(numtaps, wc);
    var hd = try allocator.alloc(f64, numtaps);
    defer allocator.free(hd);
    var i: usize = 0;
    while (i < numtaps) : (i += 1) {
        hd[i] = hd1[i] - hd2[i];
    }

    var hn = try allocator.alloc(f64, numtaps);
    const wn = try WindowFunction.kaiser(allocator, numtaps, beta);
    i = 0;
    while (i < numtaps) : (i += 1) {
        hn[i] = hd[i] * wn[i];
    }
    return hn;
}

pub fn Bandpass(numtaps: usize, cutoff: []f64, width: f64, fs: f64) ![]f64 {
    const allocator = std.heap.page_allocator;

    const atten = kaiser_atten(numtaps, width / (0.5 * fs));
    const beta = kaiser_beta(atten);
    var wc = try allocator.alloc(f64, cutoff.len);
    defer allocator.free(wc);
    var i: usize = 0;
    while (i < cutoff.len) : (i += 1) {
        cutoff[i] = cutoff[i] / (0.5 * fs);
        wc[i] = std.math.pi * cutoff[i];
    }

    const hd1 = try ideal_lp(numtaps, wc[1]);
    const hd2 = try ideal_lp(numtaps, wc[0]);
    var hd = try allocator.alloc(f64, numtaps);
    defer allocator.free(hd);
    i = 0;
    while (i < numtaps) : (i += 1) {
        hd[i] = hd1[i] - hd2[i];
    }

    var hn = try allocator.alloc(f64, numtaps);
    const wn = try WindowFunction.kaiser(allocator, numtaps, beta);
    i = 0;
    while (i < numtaps) : (i += 1) {
        hn[i] = hd[i] * wn[i];
    }
    return hn;
}

pub fn Bandstop(numtaps: usize, cutoff: []f64, width: f64, fs: f64) ![]f64 {
    const allocator = std.heap.page_allocator;

    const atten = kaiser_atten(numtaps, width / (0.5 * fs));
    const beta = kaiser_beta(atten);
    var wc = try allocator.alloc(f64, cutoff.len);
    defer allocator.free(wc);
    var i: usize = 0;
    while (i < cutoff.len) : (i += 1) {
        cutoff[i] = cutoff[i] / (0.5 * fs);
        wc[i] = std.math.pi * cutoff[i];
    }

    const hd1 = try ideal_lp(numtaps, std.math.pi);
    const hd2 = try ideal_lp(numtaps, wc[1]);
    const hd3 = try ideal_lp(numtaps, wc[0]);
    var hd = try allocator.alloc(f64, numtaps);
    defer allocator.free(hd);
    i = 0;
    while (i < numtaps) : (i += 1) {
        hd[i] = hd1[i] - (hd2[i] - hd3[i]);
    }

    var hn = try allocator.alloc(f64, numtaps);
    const wn = try WindowFunction.kaiser(allocator, numtaps, beta);
    i = 0;
    while (i < numtaps) : (i += 1) {
        hn[i] = hd[i] * wn[i];
    }
    return hn;
}

pub fn ideal_lp(numtaps: usize, wc: f64) ![]f64 {
    const allocator = std.heap.page_allocator;

    var m = try allocator.alloc(f64, numtaps);
    defer allocator.free(m);
    var i: usize = 0;
    while (i < numtaps) : (i += 1) {
        m[i] = @intToFloat(f64, i) - @intToFloat(f64, numtaps - 1) / 2 + 2.2204 * std.math.pow(f64, 10, -16);
    }

    var hd = try allocator.alloc(f64, numtaps);
    i = 0;
    while (i < numtaps) : (i += 1) {
        hd[i] = @sin(wc * m[i]) / (std.math.pi * m[i]);
    }
    return hd;
}

pub fn kaiser_atten(numtaps: usize, width: f64) f64 {
    return 2.285 * @intToFloat(f64, numtaps - 1) * std.math.pi * width + 7.95;
}

pub fn kaiser_beta(a: f64) f64 {
    var beta: f64 = undefined;
    if (a > 50) {
        beta = 0.1102 * (a - 8.7);
    } else if (a > 21) {
        beta = 0.5842 * std.math.pow(f64, (a - 21), 0.4) + 0.07886 * (a - 21);
    } else {
        beta = 0;
    }
    return beta;
}

const kaiserordResult = struct { beta: f64, numtaps: usize };

pub fn kaiserord(ripple: f64, width: f64) kaiserordResult {
    const A = @fabs(ripple);
    const beta = kaiser_beta(A);
    const numtaps = @floatToInt(usize, @ceil((A - 7.95) / 2.285 / (std.math.pi * width)) + 1);
    return kaiserordResult{ .beta = beta, .numtaps = numtaps };
}

pub fn main() !void {
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        var b = [_]f64{ 6, 7, 8, 9 };
        var ret = try Filter(&a, &b, &a);
        std.debug.print("{any}\n", .{ret});
    }
    {
        var ret = try Lowpass(5, 5, 5, 100);
        std.debug.print("{any}\n", .{ret});
    }
    {
        var ret = try Highpass(5, 5, 5, 100);
        std.debug.print("{any}\n", .{ret});
    }
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        var ret = try Bandpass(5, &a, 5, 100);
        std.debug.print("{any}\n", .{ret});
    }
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        var ret = try Bandstop(5, &a, 5, 100);
        std.debug.print("{any}\n", .{ret});
    }
    {
        var ret = try ideal_lp(5, 5);
        std.debug.print("{any}\n", .{ret});
    }
    {
        var ret = kaiser_atten(2, 2);
        std.debug.print("{}\n", .{ret});
    }
    {
        var ret = kaiser_beta(100);
        std.debug.print("{}\n", .{ret});
    }
    {
        var ret = kaiserord(100, 2);
        std.debug.print("{}\n", .{ret});
    }
}
