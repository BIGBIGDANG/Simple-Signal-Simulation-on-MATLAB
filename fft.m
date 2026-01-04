%% 实验一：FFT算法的MATLAB实现
% 功能：
% 1) 原始语音 x[n] -> 频域 X[k]
% 2) 置零 |X[k]| < Thr
% 3) 反变换重构 y[n]
% 4) FFT/IFFT 与 自编DFT/IDFT 两种实现 + 计时
% 5) 同一窗口 2x4 子图对比

clear; clc; close all;

%% ======= 初始参数 =======
N      = 1024;     % 变换点数
Thr    = 1.0;      % 频域幅值阈值
useOneSide = true; % 是否画单边幅度谱
NfftPlot  = N;     % 画频谱用的点数
% ====================================

%% 读取语音
load mtlb          % mtlb, Fs
x_all = mtlb(:);
if length(x_all) < N
    error('语音长度不足 N=%d，请减小N或换信号。', N);
end
x = x_all(1:N);
t = (0:N-1)/Fs;

%% 频率轴（用于画图）
% 画频谱用 NfftPlot 点（可零填充）
f_full = (0:NfftPlot-1) * (Fs/NfftPlot);
if useOneSide
    idx = 1:(floor(NfftPlot/2)+1);
    f   = f_full(idx);
else
    idx = 1:NfftPlot;
    f   = f_full;
end

%% ---------------- 方法A：FFT / IFFT ----------------
tic;
X_fft = fft(x, NfftPlot);                % 可零填充到 NfftPlot
X_fft_f = X_fft;
X_fft_f(abs(X_fft_f) < Thr) = 0;         % 置零小幅值频域点
y_fft_full = ifft(X_fft_f, NfftPlot);    % 重构（长度 NfftPlot）
y_fft = y_fft_full(1:N);                 % 取回前N点对齐比较
t_fft = toc;

%% ---------------- 方法B：自编DFT / IDFT（矩阵法） ----------------
% 注意：矩阵法是 O(N^2)，N 大会很慢；用于验证原理/对比即可
% 为公平起见，自编部分也按 NfftPlot 做变换，然后再取前N点
tic;
X_dft = myDFT_mat(x, NfftPlot);
X_dft_f = X_dft;
X_dft_f(abs(X_dft_f) < Thr) = 0;
y_dft_full = myIDFT_mat(X_dft_f);
y_dft = y_dft_full(1:N);
t_dft = toc;

%% ---------------- 作图（同一窗口 2x4） ----------------
figure('Name','实验一：FFT vs 自编DFT','Color','w');

% 统一谱幅坐标，避免“线条看不清”
specMax = max(abs(X_fft(idx)));
ylimSpec = [0, specMax*1.05 + eps];

% ===== 第一行：FFT/IFFT =====
subplot(2,4,1);
plot(t, x, 'LineWidth', 1); grid on;
title(sprintf('原始语音信号 x[n] (N=%d)', N));
xlabel('t/s'); ylabel('幅值');

subplot(2,4,2);
plot(f, abs(X_fft(idx)), 'LineWidth', 1); grid on;
title('FFT 幅度谱 |X|');
xlabel('f/Hz'); ylabel('|X|');
ylim(ylimSpec);

subplot(2,4,3);
plot(f, abs(X_fft_f(idx)), '.-', 'LineWidth', 1); grid on;
title(sprintf('阈值置零后 |X|<%.3g', Thr));
xlabel('f/Hz'); ylabel('|X|');
ylim(ylimSpec);

subplot(2,4,4);
plot(t, real(y_fft), 'LineWidth', 1); grid on;
title('IFFT 重构信号 y[n]');
xlabel('t/s'); ylabel('幅值');

% ===== 第二行：自编DFT/IDFT =====
subplot(2,4,5);
plot(t, x, 'LineWidth', 1); grid on;
title('原始语音信号 x[n]');
xlabel('t/s'); ylabel('幅值');

subplot(2,4,6);
plot(f, abs(X_dft(idx)), 'LineWidth', 1); grid on;
title('自编DFT 幅度谱 |X|');
xlabel('f/Hz'); ylabel('|X|');

subplot(2,4,7);
plot(f, abs(X_dft_f(idx)), '.-', 'LineWidth', 1); grid on;
title(sprintf('阈值置零后 |X|<%.3g', Thr));
xlabel('f/Hz'); ylabel('|X|');

subplot(2,4,8);
plot(t, real(y_dft), 'LineWidth', 1); grid on;
title('自编IDFT 重构信号 y[n]');
xlabel('t/s'); ylabel('幅值');

sgtitle(sprintf('运行时间：FFT/IFFT = %.6f s   |   自编DFT/IDFT = %.6f s   (N=%d, NfftPlot=%d, Thr=%.3g)', ...
    t_fft, t_dft, N, NfftPlot, Thr));

fprintf('FFT/IFFT 运行时间：%.6f s\n', t_fft);
fprintf('自编DFT/IDFT 运行时间：%.6f s\n', t_dft);

%% =============== 自编DFT/IDFT（矩阵法） ===============
function X = myDFT_mat(x, Nfft)
% DFT：X[k] = sum_{n=0}^{N-1} x[n] * exp(-j2pi kn/Nfft)
% 输入 x 允许长度<Nfft，会自动补零到Nfft
x = x(:);
N = length(x);
if N < Nfft
    x = [x; zeros(Nfft-N,1)];
elseif N > Nfft
    x = x(1:Nfft);
end

n = (0:Nfft-1).';
k = 0:Nfft-1;
W = exp(-1j*2*pi/Nfft);
WNnk = W .^ (n*k);
X = (x.' * WNnk).';   % 列向量
end

function x = myIDFT_mat(X)
% IDFT：x[n] = (1/N) * sum_{k=0}^{N-1} X[k] * exp(+j2pi kn/N)
X = X(:);
N = length(X);

n = (0:N-1).';
k = 0:N-1;
W = exp(1j*2*pi/N);
WNnk = W .^ (n*k);
x = (1/N) * (X.' * WNnk).';  % 列向量
end
