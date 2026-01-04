%% 实验二：IIR(Butterworth) 低通滤波器
clear; clc; close all;

%% ===== 设置参数 =====
Fs = 80e3;      % 采样频率(Hz)
fp = 4e3;       % 通带边界(Hz)
fsb = 20e3;     % 阻带边界(Hz)
Rp = 0.5;       % 通带最大衰减(dB)
Rs = 45;        % 阻带最小衰减(dB)
% =====================

% 归一化到奈奎斯特频率 Fs/2，范围必须在 (0,1)
Wp = fp  / (Fs/2);
Ws = fsb / (Fs/2);


if ~(Wp>0 && Wp<1 && Ws>0 && Ws<1 && Wp<Ws)
    error('频率参数不合法：需满足 0<Wp<Ws<1（注意除以Fs/2）。当前 Wp=%.3f Ws=%.3f', Wp, Ws);
end

%% 设计：最小阶数 + 截止频率
[n, Wn] = buttord(Wp, Ws, Rp, Rs);
[b, a]  = butter(n, Wn, 'low');

fprintf('n=%d, Wn=%.6f -> fc=%.2f Hz\n', n, Wn, Wn*(Fs/2));
disp('b = '); disp(b);
disp('a = '); disp(a);

%% 频响
Nfft = 4096;
[H, f] = freqz(b, a, Nfft, Fs);
magdB  = 20*log10(abs(H) + eps);
phase  = unwrap(angle(H));

%% 图
figure('Color','w','Name','IIR Butterworth LPF');
subplot(3,1,1);
stem(0:numel(b)-1, b, 'filled'); hold on; grid on;
stem(0:numel(a)-1, a, 'filled');
title('H(z) 系数'); xlabel('序号'); ylabel('系数值');
legend('b(k)','a(k)');

subplot(3,1,2);
plot(f, magdB, 'LineWidth', 1); grid on; hold on;
title('|H(f)| 幅频特性(dB)'); xlabel('f/Hz'); ylabel('dB');
xline(fp,  '--', 'fp');  xline(fsb, '--', 'fs');
yline(-Rp, ':', '-Rp');  yline(-Rs, ':', '-Rs');
[~, ip] = min(abs(f-fp));  [~, is] = min(abs(f-fsb));
plot(f(ip), magdB(ip), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(f(is), magdB(is), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);

subplot(3,1,3);
plot(f, phase, 'LineWidth', 1); grid on; hold on;
title('相频特性(展开)'); xlabel('f/Hz'); ylabel('rad');
xline(fp,  '--', 'fp');  xline(fsb, '--', 'fs');
