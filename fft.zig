//!为快速傅立叶变换和逆变换提供方法

const std = @import("std");

///快速傅立叶变换
pub fn fft(real: []f64, imag: []f64, truncate: bool) ![2][]f64 {
    const allocator = std.heap.page_allocator;

    const data_length = real.len;
    if (data_length <= 1) return [2][]f64{ try allocator.dupe(f64, real), try allocator.dupe(f64, imag) };

    const n = if (truncate) previous_pow2(@intCast(isize, data_length)) else next_pow2(@intCast(isize, data_length));

    var doubleReal = try allocator.alloc(f64, @intCast(usize, n));
    defer allocator.free(doubleReal);
    var doubleImag = try allocator.alloc(f64, @intCast(usize, n));
    defer allocator.free(doubleImag);

    std.mem.copy(f64, doubleReal, real);
    std.mem.copy(f64, doubleImag, imag);

    try pow2_fft(doubleReal, doubleImag);

    var retReal = try allocator.alloc(f64, @intCast(usize, n));
    var retImag = try allocator.alloc(f64, @intCast(usize, n));

    var i: usize = 0;
    while (i < n) : (i += 1) {
        retReal[i] = doubleReal[i];
        retImag[i] = doubleImag[i];
    }
    return [2][]f64{ retReal, retImag };
}

///快速傅立叶逆变换
pub fn ifft(real: []f64, imag: []f64, truncate: bool) ![2][]f64 {
    const allocator = std.heap.page_allocator;

    const data_length = real.len;
    if (data_length <= 1) return [2][]f64{ try allocator.dupe(f64, real), try allocator.dupe(f64, imag) };

    const n = if (truncate) previous_pow2(@intCast(isize, data_length)) else next_pow2(@intCast(isize, data_length));

    var doubleReal = try allocator.alloc(f64, @intCast(usize, n));
    defer allocator.free(doubleReal);
    var doubleImag = try allocator.alloc(f64, @intCast(usize, n));
    defer allocator.free(doubleImag);

    std.mem.copy(f64, doubleReal, real);
    std.mem.copy(f64, doubleImag, imag);

    try pow2_ifft(doubleReal, doubleImag);

    var retReal = try allocator.alloc(f64, @intCast(usize, n));
    var retImag = try allocator.alloc(f64, @intCast(usize, n));

    var i: usize = 0;
    while (i < n) : (i += 1) {
        retReal[i] = doubleReal[i];
        retImag[i] = doubleImag[i];
    }
    return [2][]f64{ retReal, retImag };
}

///快速傅立叶变换
///输入数据长度必须是2的n次幂, 且n >= 1
pub fn pow2_fft(real: []f64, imag: []f64) !void {
    const n = real.len;
    bitrp(real, imag);
    //计算 1 的前 n / 2 个 n 次方根的共轭复数 W'j = wr [j] + i * wi [j] , j = 0, 1, ... , n / 2 - 1
    const half_n = (n >> 1);

    const allocator = std.heap.page_allocator;

    var wr = try allocator.alloc(f64, half_n);
    defer allocator.free(wr);
    var wi = try allocator.alloc(f64, half_n);
    defer allocator.free(wi);

    const arg = -2 * std.math.pi / @intToFloat(f64, n);
    const cos = @cos(arg);
    const sin = @sin(arg);
    wr[0] = 1;
    wi[0] = 0;
    var i: usize = 1;
    while (i < half_n) : (i += 1) {
        wr[i] = wr[i - 1] * cos - wi[i - 1] * sin;
        wi[i] = wr[i - 1] * sin + wi[i - 1] * cos;
    }
    var m: usize = 2;
    while (m <= n) : (m <<= 1) {
        var half_m = (m >> 1);
        var j: usize = 0;
        while (j < n) : (j += m) {
            i = 0;
            while (i < half_m) : (i += 1) {
                var index1 = i + j;
                var index2 = index1 + half_m;
                var t = n / m * i;
                if (index1 >= n or index2 >= n) continue; //不加此代码报错
                var tr = wr[t] * real[index2] - wi[t] * imag[index2];
                var ti = wr[t] * imag[index2] + wi[t] * real[index2];
                var ur = real[index1];
                var ui = imag[index1];
                real[index1] = ur + tr;
                imag[index1] = ui + ti;
                real[index2] = ur - tr;
                imag[index2] = ui - ti;
            }
        }
    }
}

///快速傅立叶逆变换
///输入数据长度必须是2的n次幂, 且n >= 1
pub fn pow2_ifft(real: []f64, imag: []f64) !void {
    const n = real.len;
    bitrp(real, imag);
    //计算 1 的前 n / 2 个 n 次方根的共轭复数 W'j = wr [j] + i * wi [j] , j = 0, 1, ... , n / 2 - 1
    const half_n = (n >> 1);

    const allocator = std.heap.page_allocator;

    var wr = try allocator.alloc(f64, half_n);
    defer allocator.free(wr);
    var wi = try allocator.alloc(f64, half_n);
    defer allocator.free(wi);

    const arg = 2 * std.math.pi / @intToFloat(f64, n);
    const cos = @cos(arg);
    const sin = @sin(arg);
    wr[0] = 1;
    wi[0] = 0;
    var i: usize = 1;
    while (i < half_n) : (i += 1) {
        wr[i] = wr[i - 1] * cos - wi[i - 1] * sin;
        wi[i] = wr[i - 1] * sin + wi[i - 1] * cos;
    }
    var m: usize = 2;
    while (m <= n) : (m <<= 1) {
        var half_m = (m >> 1);
        var j: usize = 0;
        while (j < n) : (j += m) {
            i = 0;
            while (i < half_m) : (i += 1) {
                var index1 = i + j;
                var index2 = index1 + half_m;
                var t = n / m * i;
                if (index1 >= n or index2 >= n) continue; //不加此代码报错
                var tr = wr[t] * real[index2] - wi[t] * imag[index2];
                var ti = wr[t] * imag[index2] + wi[t] * real[index2];
                var ur = real[index1];
                var ui = imag[index1];
                real[index1] = ur + tr;
                imag[index1] = ui + ti;
                real[index2] = ur - tr;
                imag[index2] = ui - ti;
            }
        }
    }
    i = 0;
    while (i < n) : (i += 1) {
        real[i] /= @intToFloat(f64, n);
        imag[i] /= @intToFloat(f64, n);
    }
}

///位反转置换(Bit-reversal Permutation)
fn bitrp(real: []f64, imag: []f64) void {
    const n = real.len;
    var p: usize = 0;
    var k: usize = 1;
    while (k < n) {
        k <<= 1;
        p += 1;
    }

    var i: usize = 0;
    while (i < n) : (i += 1) {
        var a = i;
        var b: usize = 0;
        var j: usize = 0;
        while (j < p) : (j += 1) {
            b = (b << 1) + (a & 0x1); // b = b * 2 + a % 2;
            a >>= 1; // a = a / 2;
        }
        if (b >= n) continue; //不加此代码报错
        if (b > i) {
            var t = real[i];
            real[i] = real[b];
            real[b] = t;

            t = imag[i];
            imag[i] = imag[b];
            imag[b] = t;
        }
    }
}

///返回大于等于n的最小2次幂(最大1G, 即0x40000000)
pub fn next_pow2(n: isize) isize {
    var ret: usize = 1;
    while (ret < n) {
        ret <<= 1;
    }
    if (ret <= 0) return 0x40000000;
    return @intCast(isize, ret);
}

///返回小于等于n的最大2次幂
pub fn previous_pow2(n: isize) isize {
    var ret: usize = 1;
    while ((ret << 1) <= n) {
        ret <<= 1;
    }
    return @intCast(isize, ret);
}

pub fn main() !void {
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        // var a = [_]f64{1};
        var b = [_]f64{ 6, 7, 8, 9 };
        // var b = [_]f64{6};
        var ret = try fft(&a, &b, true);
        std.debug.print("{any}\n{any}\n", .{ a, b });
        std.debug.print("{any}\n{any}\n", .{ ret[0], ret[1] });
    }
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        // var a = [_]f64{1};
        var b = [_]f64{ 6, 7, 8, 9 };
        // var b = [_]f64{6};
        var ret = try ifft(&a, &b, true);
        std.debug.print("{any}\n{any}\n", .{ a, b });
        std.debug.print("{any}\n{any}\n", .{ ret[0], ret[1] });
    }
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        var b = [_]f64{ 6, 7, 8, 9 };
        try pow2_fft(&a, &b);
        std.debug.print("{any}\n{any}\n", .{ a, b });
    }
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        var b = [_]f64{ 6, 7, 8, 9 };
        try pow2_ifft(&a, &b);
        std.debug.print("{any}\n{any}\n", .{ a, b });
    }
    {
        var a = [_]f64{ 1, 2, 3, 4, 5 };
        var b = [_]f64{ 6, 7, 8, 9, 10 };
        bitrp(&a, &b);
        std.debug.print("{any}\n{any}\n", .{ a, b });
    }
    {
        var ret = next_pow2(100);
        std.debug.print("{}\n", .{ret});
    }
    {
        var ret = previous_pow2(100);
        std.debug.print("{}\n", .{ret});
    }
}
