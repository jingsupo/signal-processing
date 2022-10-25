//!提供各种窗函数

const std = @import("std");

///平顶窗
pub fn flattop(allocator: std.mem.Allocator, n: usize, periodic: bool) ![]f64 {
    var w = std.ArrayList(f64).init(allocator);
    defer w.deinit();

    if (n <= 0) return w.toOwnedSlice();
    if (n == 1) {
        try w.append(1.0);
        return w.toOwnedSlice();
    }

    var para = if (periodic) (2.0 * std.math.pi / @intToFloat(f64, n)) else (2.0 * std.math.pi / @intToFloat(f64, n - 1));
    const a0 = 0.21557895;
    const a1 = 0.41663158;
    const a2 = 0.277263158;
    const a3 = 0.083578947;
    const a4 = 0.006947368;
    var i: usize = 0;
    while (i < n) : (i += 1) {
        try w.append(a0 -
            a1 * @cos(@intToFloat(f64, i) * para) +
            a2 * @cos(@intToFloat(f64, 2 * i) * para) -
            a3 * @cos(@intToFloat(f64, 3 * i) * para) +
            a4 * @cos(@intToFloat(f64, 4 * i) * para));
    }
    return w.toOwnedSlice();
}

///海明窗
pub fn hamming(allocator: std.mem.Allocator, n: usize, periodic: bool) ![]f64 {
    var w = std.ArrayList(f64).init(allocator);
    defer w.deinit();

    if (n <= 0) return w.toOwnedSlice();
    if (n == 1) {
        try w.append(1.0);
        return w.toOwnedSlice();
    }

    var para = if (periodic) (2.0 * std.math.pi / @intToFloat(f64, n)) else (2.0 * std.math.pi / @intToFloat(f64, n - 1));
    var i: usize = 0;
    while (i < n) : (i += 1) {
        try w.append(0.54 - 0.46 * @cos(@intToFloat(f64, i) * para));
    }
    return w.toOwnedSlice();
}

///汉宁窗
pub fn hann(allocator: std.mem.Allocator, n: usize, periodic: bool) ![]f64 {
    var w = std.ArrayList(f64).init(allocator);
    defer w.deinit();

    if (n <= 0) return w.toOwnedSlice();
    if (n == 1) {
        try w.append(1.0);
        return w.toOwnedSlice();
    }

    var para = if (periodic) (2.0 * std.math.pi / @intToFloat(f64, n)) else (2.0 * std.math.pi / @intToFloat(f64, n - 1));
    var i: usize = 0;
    while (i < n) : (i += 1) {
        try w.append(0.5 - 0.5 * @cos(@intToFloat(f64, i) * para));
    }
    return w.toOwnedSlice();
}

///布莱克曼窗
pub fn blackman(allocator: std.mem.Allocator, n: usize, periodic: bool) ![]f64 {
    var w = std.ArrayList(f64).init(allocator);
    defer w.deinit();

    if (n <= 0) return w.toOwnedSlice();
    if (n == 1) {
        try w.append(1.0);
        return w.toOwnedSlice();
    }

    var para = if (periodic) (2.0 * std.math.pi / @intToFloat(f64, n)) else (2.0 * std.math.pi / @intToFloat(f64, n - 1));
    var i: usize = 0;
    while (i < n) : (i += 1) {
        try w.append(0.42 -
            0.5 * @cos(@intToFloat(f64, i) * para) +
            0.08 * @cos(@intToFloat(f64, 2 * i) * para));
    }
    return w.toOwnedSlice();
}

///巴特莱特窗
pub fn bartlett(allocator: std.mem.Allocator, n: usize) ![]f64 {
    var w = std.ArrayList(f64).init(allocator);
    defer w.deinit();

    if (n <= 0) return w.toOwnedSlice();
    if (n == 1) {
        try w.append(1.0);
        return w.toOwnedSlice();
    }

    const half = (n - 1) >> 1;
    {
        var i: usize = 0;
        while (i <= half) : (i += 1) {
            try w.append(@intToFloat(f64, i * 2) / @intToFloat(f64, n - 1));
        }
    }
    {
        var i: usize = half + 1;
        while (i < n) : (i += 1) {
            try w.append(2 - @intToFloat(f64, i * 2) / @intToFloat(f64, n - 1));
        }
    }
    return w.toOwnedSlice();
}

///三角窗(费杰窗)
pub fn triang(allocator: std.mem.Allocator, n: usize) ![]f64 {
    var w = std.ArrayList(f64).init(allocator);
    defer w.deinit();

    if (n <= 0) return w.toOwnedSlice();

    const odd = (n & 1) == 1;
    const half = if (odd) (n + 1) / 2 else n / 2;
    if (odd) {
        {
            var i: usize = 0;
            while (i < half) : (i += 1) {
                try w.append(@intToFloat(f64, 2 * (i + 1)) / @intToFloat(f64, n + 1));
            }
        }
        {
            var i: usize = half;
            while (i < n) : (i += 1) {
                try w.append(2 - @intToFloat(f64, 2 * (i + 1)) / @intToFloat(f64, n + 1));
            }
        }
    } else {
        {
            var i: usize = 0;
            while (i < half) : (i += 1) {
                try w.append(@intToFloat(f64, (2 * (i + 1) - 1)) / @intToFloat(f64, n));
            }
        }
        {
            var i: usize = half;
            while (i < n) : (i += 1) {
                try w.append(2 - @intToFloat(f64, (2 * (i + 1) - 1)) / @intToFloat(f64, n));
            }
        }
    }
    return w.toOwnedSlice();
}

///高斯窗
pub fn gaussian(allocator: std.mem.Allocator, n: usize, std_: f64) ![]f64 {
    var w = std.ArrayList(f64).init(allocator);
    defer w.deinit();

    if (n <= 0) return w.toOwnedSlice();
    if (n == 1) {
        try w.append(1.0);
        return w.toOwnedSlice();
    }

    const half = 0.5 * @intToFloat(f64, n - 1);
    const sig2 = std_ * std_;
    var i: usize = 0;
    while (i < n) : (i += 1) {
        try w.append(@exp(-0.5 * ((@intToFloat(f64, i) - half) * (@intToFloat(f64, i) - half) / sig2)));
    }
    return w.toOwnedSlice();
}

///凯撒窗
pub fn kaiser(allocator: std.mem.Allocator, n: usize, beta: f64) ![]f64 {
    var w = std.ArrayList(f64).init(allocator);
    defer w.deinit();

    const bes = first_modified_Bessel(0, beta);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        const x = 1.0 - 2.0 * @intToFloat(f64, i) / @intToFloat(f64, n - 1);
        const para = beta * @sqrt(1 - x * x);
        try w.append(first_modified_Bessel(0, para) / bes);
    }
    return w.toOwnedSlice();
}

///第一类修正贝塞尔函数
fn first_modified_Bessel(n_: isize, x: f64) f64 {
    var n = n_;
    if (n < 0) n = -n;

    var p: f64 = 0;
    var q: f64 = 0;
    const a = [7]f64{ 1.0, 3.5156229, 3.0899424, 1.2067492, 0.2659732, 0.0360768, 0.0045813 };
    const b = [7]f64{ 0.5, 0.87890594, 0.51498869, 0.15084934, 0.02658773, 0.00301532, 0.00032411 };
    const c = [9]f64{ 0.39894228, 0.01328592, 0.00225319, -0.00157565, 0.00916281, -0.02057706, 0.02635537, -0.01647633, 0.00392377 };
    const d = [9]f64{ 0.39894228, -0.03988024, -0.00362018, 0.00163801, -0.01031555, 0.02282967, -0.02895312, 0.01787654, -0.00420059 };

    const t = @fabs(x);
    if (n != 1) {
        if (t < 3.75) {
            const y = (x / 3.75) * (x / 3.75);
            p = a[6];
            var i: isize = 5;
            while (i >= 0) : (i -= 1) {
                p = p * y + a[@intCast(usize, i)];
            }
        } else {
            const y = 3.75 / t;
            p = c[8];
            var i: isize = 7;
            while (i >= 0) : (i -= 1) {
                p = p * y + c[@intCast(usize, i)];
            }
            p *= @exp(t) / @sqrt(t);
        }
    }

    if (n == 0) return p;

    q = p;
    if (t < 3.75) {
        const y = (x / 3.75) * (x / 3.75);
        p = b[6];
        var i: isize = 5;
        while (i >= 0) : (i -= 1) {
            p = p * y + b[@intCast(usize, i)];
        }
        p *= t;
    } else {
        const y = 3.75 / t;
        p = d[8];
        var i: isize = 7;
        while (i >= 0) : (i -= 1) {
            p = p * y + d[@intCast(usize, i)];
        }
        p *= @exp(t) / @sqrt(t);
    }

    if (x < 0) p = -p;

    if (n == 1) return p;

    if (x == 0) return 0;

    const z = 2.0 / t;
    var bt: f64 = 0;
    var b1: f64 = 1;
    var b0: f64 = 0;
    const m: usize = @intCast(usize, 2 * (n + @floatToInt(isize, @sqrt(40.0 * @intToFloat(f64, n)))));

    const powPlus = std.math.pow(f64, 10, 10);
    const powMinus = std.math.pow(f64, 10, -10);

    var i: usize = m;
    while (i > 0) : (i -= 1) {
        p = b0 + @intToFloat(f64, i) * z * b1;
        b0 = b1;
        b1 = p;
        if (@fabs(b1) > powPlus) {
            bt *= powMinus;
            b0 *= powMinus;
            b1 *= powMinus;
        }
        if (i == n) bt = b0;
    }
    p = bt * q / b1;
    if ((x < 0) and ((n & 1) == 1)) p = -p;
    return p;
}

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    {
        var w = try flattop(allocator, 5, true);
        var sum: f64 = 0;
        for (w) |v| {
            std.debug.print("\n{}\n", .{v});
            sum += v;
        }
        std.debug.print("\n{}\n", .{sum});
    }
    {
        var w = try hamming(allocator, 5, false);
        var sum: f64 = 0;
        for (w) |v| {
            std.debug.print("\n{}\n", .{v});
            sum += v;
        }
        std.debug.print("\n{}\n", .{sum});
    }
    {
        var w = try hann(allocator, 5, false);
        var sum: f64 = 0;
        for (w) |v| {
            std.debug.print("\n{}\n", .{v});
            sum += v;
        }
        std.debug.print("\n{}\n", .{sum});
    }
    {
        var w = try blackman(allocator, 5, false);
        var sum: f64 = 0;
        for (w) |v| {
            std.debug.print("\n{}\n", .{v});
            sum += v;
        }
        std.debug.print("\n{}\n", .{sum});
    }
    {
        var w = try bartlett(allocator, 5);
        var sum: f64 = 0;
        for (w) |v| {
            std.debug.print("\n{}\n", .{v});
            sum += v;
        }
        std.debug.print("\n{}\n", .{sum});
    }
    {
        var w = try triang(allocator, 5);
        var sum: f64 = 0;
        for (w) |v| {
            std.debug.print("\n{}\n", .{v});
            sum += v;
        }
        std.debug.print("\n{}\n", .{sum});
    }
    {
        var w = try gaussian(allocator, 5, 7);
        var sum: f64 = 0;
        for (w) |v| {
            std.debug.print("\n{}\n", .{v});
            sum += v;
        }
        std.debug.print("\n{}\n", .{sum});
    }
    {
        var w = try kaiser(allocator, 5, 5);
        var sum: f64 = 0;
        for (w) |v| {
            std.debug.print("\n{}\n", .{v});
            sum += v;
        }
        std.debug.print("\n{}\n", .{sum});
    }
}
