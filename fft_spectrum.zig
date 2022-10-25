//!为画fft幅值谱和功率谱提供方法

const std = @import("std");

const FFT = @import("fft.zig");

///线性平均法计算fft幅值谱图
///<param name="data">输入数据</param>
///<param name="frequency">采样频率</param>
///<param name="line">谱线数</param>
///<param name="averageCount">平均次数(即块数)</param>
///<param name="windowFunction">窗函数</param>
///<returns></returns>
pub fn AmplitudeSpectrumLinearAverage(data: []f64, frequency: f64, line: usize, averageCount: usize, windowFunction: ?[]f64) ![2][]f64 {
    const allocator = std.heap.page_allocator;

    const amplitudes = try fft_amplitude(data, line, averageCount, windowFunction);
    const n = amplitudes.len / 2 + 1;
    var x = try allocator.alloc(f64, n);
    var y = try allocator.alloc(f64, n);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        x[i] = @intToFloat(f64, i) * (frequency / @intToFloat(f64, amplitudes.len));
    }
    std.mem.copy(f64, y, amplitudes);
    return [2][]f64{ x, y };
}

///线性平均法计算fft功率谱图
///<param name="data">输入数据</param>
///<param name="frequency">采样频率</param>
///<param name="line">谱线数</param>
///<param name="averageCount">平均次数(即块数)</param>
///<param name="windowFunction">窗函数</param>
///<returns></returns>
pub fn PowerSpectrumLinearAverage(data: []f64, frequency: f64, line: usize, averageCount: usize, windowFunction: ?[]f64) ![2][]f64 {
    const allocator = std.heap.page_allocator;

    const powers = try fft_power(data, line, averageCount, windowFunction);
    const n = powers.len / 2 + 1;
    var x = try allocator.alloc(f64, n);
    var y = try allocator.alloc(f64, n);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        x[i] = @intToFloat(f64, i) * frequency;
    }
    std.mem.copy(f64, y, powers);
    return [2][]f64{ x, y };
}

///fft幅值谱图
///<param name="data">输入数据</param>
///<param name="frequency">采样频率</param>
///<returns></returns>
pub fn AmplitudeSpectrum(data: []f64, frequency: f64) ![2][]f64 {
    const allocator = std.heap.page_allocator;

    var tmp = try allocator.alloc(f64, data.len);
    defer allocator.free(tmp);
    const fft_result = try FFT.fft(data, tmp, true);
    const real = fft_result[0];
    const imag = fft_result[1];
    const n = real.len / 2 + 1;
    var x = try allocator.alloc(f64, n);
    var y = try allocator.alloc(f64, n);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        x[i] = @intToFloat(f64, i) * (frequency / @intToFloat(f64, real.len));
        y[i] = 2 * @sqrt(real[i] * real[i] + imag[i] * imag[i]) / @intToFloat(f64, data.len); //matlab帮助里的画法
        // y[i] = 2 * @sqrt(real[i] * real[i] + imag[i] * imag[i]); //论坛中的画法
    }
    return [2][]f64{ x, y };
}

///fft功率谱图
///<param name="data">输入数据</param>
///<param name="frequency">采样频率</param>
///<returns></returns>
pub fn PowerSpectrum(data: []f64, frequency: f64) ![2][]f64 {
    const allocator = std.heap.page_allocator;

    var tmp = try allocator.alloc(f64, data.len);
    defer allocator.free(tmp);
    const fft_result = try FFT.fft(data, tmp, true);
    const real = fft_result[0];
    const imag = fft_result[1];
    const n = real.len / 2 + 1;
    var x = try allocator.alloc(f64, n);
    var y = try allocator.alloc(f64, n);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        x[i] = @intToFloat(f64, i) * frequency;
        // y[i] = 2 * @sqrt(real[i] * real[i] + imag[i] * imag[i]) / @intToFloat(f64, data.len); //matlab帮助里的画法
        y[i] = real[i] * real[i] + imag[i] * imag[i]; //论坛中的画法
    }
    return [2][]f64{ x, y };
}

///线性平均法计算fft幅值(分块, 加窗, 每块做fft, 每块中每个点计算幅值, 按照块数计算相应点的幅值平均值)
///<param name="data">输入数据</param>
///<param name="line">谱线数</param>
///<param name="averageCount">平均次数(即块数)</param>
///<param name="windowFunction">窗函数</param>
///<returns>幅值</returns>
fn fft_amplitude(data: []f64, line: usize, averageCount: usize, windowFunction: ?[]f64) ![]f64 {
    const allocator = std.heap.page_allocator;

    const blockSize = @floatToInt(usize, @intToFloat(f64, line) * 2.56);
    const blockCount = averageCount;
    const total = blockSize * blockCount;

    var innerData = try allocator.alloc(f64, total);
    defer allocator.free(innerData);
    std.mem.copy(f64, innerData, data);

    var ret = try allocator.alloc(f64, blockSize);
    var i: usize = 0;
    while (i < blockCount) : (i += 1) {
        var real = try allocator.alloc(f64, blockSize);
        defer allocator.free(real);
        var imag = try allocator.alloc(f64, blockSize);
        defer allocator.free(imag);
        const sourceIndex = i * blockSize;
        for (innerData[sourceIndex..(sourceIndex + blockSize)]) |value, index|
            real[index] = value;
        if (windowFunction != null) {
            var k: usize = 0;
            while (k < blockSize) : (k += 1) {
                real[k] *= windowFunction.?[k];
            }
        }

        try FFT.pow2_fft(real, imag);

        var j: usize = 0;
        while (j < blockSize) : (j += 1) {
            var amplitude = 2 * @sqrt(real[j] * real[j] + imag[j] * imag[j]) / @intToFloat(f64, blockSize);
            ret[j] += amplitude / @intToFloat(f64, blockCount);
        }
    }
    return ret;
}

///线性平均法计算fft功率(分块, 加窗, 每块做fft, 每块中每个点计算幅值, 按照块数计算相应点的功率平均值)
///<param name="data">输入数据</param>
///<param name="line">谱线数</param>
///<param name="averageCount">平均次数(即块数)</param>
///<param name="windowFunction">窗函数</param>
///<returns>功率</returns>
fn fft_power(data: []f64, line: usize, averageCount: usize, windowFunction: ?[]f64) ![]f64 {
    const allocator = std.heap.page_allocator;

    const blockSize = @floatToInt(usize, @intToFloat(f64, line) * 2.56);
    const blockCount = averageCount;
    const total = blockSize * blockCount;

    var innerData = try allocator.alloc(f64, total);
    defer allocator.free(innerData);
    std.mem.copy(f64, innerData, data);

    var ret = try allocator.alloc(f64, blockSize);
    var i: usize = 0;
    while (i < blockCount) : (i += 1) {
        var real = try allocator.alloc(f64, blockSize);
        defer allocator.free(real);
        var imag = try allocator.alloc(f64, blockSize);
        defer allocator.free(imag);
        const sourceIndex = i * blockSize;
        for (innerData[sourceIndex..(sourceIndex + blockSize)]) |value, index|
            real[index] = value;
        if (windowFunction != null) {
            var k: usize = 0;
            while (k < blockSize) : (k += 1) {
                real[k] *= windowFunction.?[k];
            }
        }

        try FFT.pow2_fft(real, imag);

        var j: usize = 0;
        while (j < blockSize) : (j += 1) {
            var power = real[j] * real[j] + imag[j] * imag[j];
            ret[j] += power / @intToFloat(f64, blockCount);
        }
    }
    return ret;
}

pub fn main() !void {
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        var b = [_]f64{ 6, 7, 8, 9 };
        var ret = try AmplitudeSpectrumLinearAverage(&a, 100, 1, 2, &b);
        std.debug.print("{any}\n", .{ret});
    }
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        var b = [_]f64{ 6, 7, 8, 9 };
        var ret = try PowerSpectrumLinearAverage(&a, 100, 1, 2, &b);
        std.debug.print("{any}\n", .{ret});
    }
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        var ret = try AmplitudeSpectrum(&a, 100);
        std.debug.print("{any}\n", .{ret});
    }
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        var ret = try PowerSpectrum(&a, 100);
        std.debug.print("{any}\n", .{ret});
    }
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        var b = [_]f64{ 6, 7, 8, 9 };
        var ret = fft_amplitude(&a, 1, 2, &b);
        std.debug.print("{any}\n", .{ret});
    }
    {
        var a = [_]f64{ 1, 2, 3, 4 };
        var b = [_]f64{ 6, 7, 8, 9 };
        var ret = fft_power(&a, 1, 2, &b);
        std.debug.print("{any}\n", .{ret});
    }
}
