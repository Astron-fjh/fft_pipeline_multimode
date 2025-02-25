%% Description:
%   本代码实现了512点FFT和IFFT的定点化处理，并与标准FFT和IFFT进行对比。
%   首先，生成一个时域信号并进行FFT，之后将信号转换为16位定点表示并保存为文本文件，用作RTL的输入。
%   随后，使用自定义的定点化FFT算法进行计算，并保存结果。同时，计算标准IFFT并进行定点化处理。
%   最后，代码可视化标准FFT、定点化FFT、标准IFFT和定点化IFFT的结果，并与RTL输出进行对比。
%
%   This code implements the fixed-point processing of 512-point FFT and IFFT, and compares the results with standard FFT and IFFT.
%   First, a time-domain signal is generated and processed with FFT. The signal is then converted to a 16-bit fixed-point representation and saved as a text file to be used as input for RTL.
%   Next, the custom fixed-point FFT algorithm is used for computation, and the results are saved. The standard IFFT is also computed and processed with fixed-point representation.
%   Finally, the code visualizes the results of the standard FFT, fixed-point FFT, standard IFFT, and fixed-point IFFT, and compares them with the RTL output.

%% Author(s):
%   Astron-fjh

% 初始化
clear all
clc
N = 512;    % 512点
t = 1:1:N;
X_re = cos(1/3*pi*t);
X_im = zeros(1, N);
X = complex(X_re, X_im);
X_fft = (fft(X));

% 时域信号转换为二进制补码并保存
X_re_16bit = floor(X_re * (2^7));   % 定点化为16位整数，放大至 2^7
X_im_16bit = floor(X_im * (2^7));

% 保证转换为有符号16位整数
X_re_int16 = int16(X_re_16bit);     % 转换为有符号16位整数
X_im_int16 = int16(X_im_16bit);

% 打开文件并写入数据
fd_re = fopen("./512点/x_re512.txt", 'wb');
fd_im = fopen("./512点/x_im512.txt", 'wb');
for i = 1:N
    % 将16位有符号整数转为补码并保存
    fprintf(fd_re, '%04x\r\n', typecast(X_re_int16(i), 'uint16'));
    fprintf(fd_im, '%04x\r\n', typecast(X_im_int16(i), 'uint16'));
end
fclose(fd_re);
fclose(fd_im);

% 频域信号（FFT结果）转换为二进制补码并保存
X_fft_re_16bit = floor(real(X_fft) * (2^7));  % 放大至2^7
X_fft_im_16bit = floor(imag(X_fft) * (2^7));

% 转换为16位有符号整数
X_fft_re_int16 = int16(X_fft_re_16bit);
X_fft_im_int16 = int16(X_fft_im_16bit);

fd_re = fopen("./512点/x_fft512_re.txt", 'wb');
fd_im = fopen("./512点/x_fft512_im.txt", 'wb');
for i = 1:N
    fprintf(fd_re, '%04x\r\n', typecast(X_fft_re_int16(i), 'uint16'));
    fprintf(fd_im, '%04x\r\n', typecast(X_fft_im_int16(i), 'uint16'));
end
fclose(fd_re);
fclose(fd_im);

% 定点化FFT计算
X_16bit = complex(X_re_16bit, X_im_16bit);
X_sim_fft = fft_radix2(X_16bit, N);

% 转换并保存定点化FFT计算结果
X_sim_fft_re_16bit = floor(real(X_sim_fft));
X_sim_fft_im_16bit = floor(imag(X_sim_fft));

X_sim_fft_re_int16 = int16(X_sim_fft_re_16bit);
X_sim_fft_im_int16 = int16(X_sim_fft_im_16bit);

fd_re = fopen("./512点/x_sim_fft512_re.txt", 'wb');
fd_im = fopen("./512点/x_sim_fft512_im.txt", 'wb');
for i = 1:N
    fprintf(fd_re, '%04x\r\n', typecast(X_sim_fft_re_int16(i), 'uint16'));
    fprintf(fd_im, '%04x\r\n', typecast(X_sim_fft_im_int16(i), 'uint16'));
end
fclose(fd_re);
fclose(fd_im);

X_sim_fft_re = floor(real(X_sim_fft)) / (2^7);
X_sim_fft_im = floor(imag(X_sim_fft)) / (2^7);
X_sim_fft = complex(X_sim_fft_re, X_sim_fft_im);

% 可视化比较
figure;
set(gcf, 'Position', [100, 100, 1200, 900]);
subplot(2, 3, 1); 
plot(abs(X_fft), 'LineWidth', 1.5);
title('标准512点FFT结果（幅度）');
xlabel('频率索引');
ylabel('幅度');
grid on;

subplot(2, 3, 2); 
plot(abs(X_sim_fft), 'LineWidth', 1.5);
title('定点化512点FFT结果（幅度）');
xlabel('频率索引');
ylabel('幅度');
grid on;

%%
% 读取RTL的输出
re_file = './512点/y_fft512_re.txt';
im_file = './512点/y_fft512_im.txt';

% 读取实部和虚部数据
y_fft512_re = [];
y_fft512_im = [];

% 读取实部数据
fid_re = fopen(re_file, 'r');
tline = fgets(fid_re);
while ischar(tline)
    y_fft512_re = [y_fft512_re; typecast(uint16(hex2dec(tline(1:4))), 'int16')];
    tline = fgets(fid_re);
end
fclose(fid_re);

% 读取虚部数据
fid_im = fopen(im_file, 'r');
tline = fgets(fid_im);
while ischar(tline)
    y_fft512_im = [y_fft512_im; typecast(uint16(hex2dec(tline(1:4))), 'int16')];
    tline = fgets(fid_im);
end
fclose(fid_im);

% 合并实部和虚部为复数数组
y_fft = complex(double(y_fft512_re) / (2^7), double(y_fft512_im) / (2^7));

subplot(2, 3, 3); 
plot(abs(y_fft), 'LineWidth', 1.5);
title('定点化RTL 512点FFT结果（幅度）');
xlabel('频率索引');
ylabel('幅度');
grid on;

%%
% 标准IFFT计算
X_ifft = ifft(X);

% 标准IFFT的定点化处理（与FFT相同）
X_ifft_re_16bit = floor(real(X_ifft) * (2^7));  % 放大至2^7
X_ifft_im_16bit = floor(imag(X_ifft) * (2^7));  % 放大至2^7

% 转换为16位有符号整数
X_ifft_re_int16 = int16(X_ifft_re_16bit);
X_ifft_im_int16 = int16(X_ifft_im_16bit);

% 保存标准IFFT的定点化结果
fd_re = fopen("./512点/x_ifft512_re.txt", 'wb');
fd_im = fopen("./512点/x_ifft512_im.txt", 'wb');
for i = 1:N
    fprintf(fd_re, '%04x\r\n', typecast(X_ifft_re_int16(i), 'uint16'));
    fprintf(fd_im, '%04x\r\n', typecast(X_ifft_im_int16(i), 'uint16'));
end
fclose(fd_re);
fclose(fd_im);

% 定点化IFFT计算
% 对于IFFT，先将FFT的输入数据取共轭，再对FFT结果取共轭，并在蝶形运算后进行右移6位
X_16bit_ifft = complex(X_re_16bit, X_im_16bit);
X_16bit_ifft = conj(X_16bit_ifft);  % 取共轭
X_sim_ifft = conj(fft_radix2(X_16bit_ifft, N)) / (2^9); 

% 保存定点化IFFT结果
X_sim_ifft_re_16bit = floor(real(X_sim_ifft));
X_sim_ifft_im_16bit = floor(imag(X_sim_ifft));

X_sim_ifft_re_int16 = int16(X_sim_ifft_re_16bit);
X_sim_ifft_im_int16 = int16(X_sim_ifft_im_16bit);

fd_re = fopen("./512点/x_sim_ifft512_re.txt", 'wb');
fd_im = fopen("./512点/x_sim_ifft512_im.txt", 'wb');
for i = 1:N
    fprintf(fd_re, '%04x\r\n', typecast(X_sim_ifft_re_int16(i), 'uint16'));
    fprintf(fd_im, '%04x\r\n', typecast(X_sim_ifft_im_int16(i), 'uint16'));
end
fclose(fd_re);
fclose(fd_im);

X_sim_ifft = floor(X_sim_ifft) / (2^7);

% 可视化标准IFFT与定点化IFFT对比
subplot(2, 3, 4); 
plot(abs(X_ifft), 'LineWidth', 1.5);
title('标准512点IFFT结果（幅度）');
xlabel('频率索引');
ylabel('幅度');
grid on;

subplot(2, 3, 5); 
plot(abs(X_sim_ifft), 'LineWidth', 1.5);
title('定点化512点IFFT结果（幅度）');
xlabel('频率索引');
ylabel('幅度');
grid on;

%%
% 读取RTL的输出
re_file = './512点/y_ifft512_re.txt';
im_file = './512点/y_ifft512_im.txt';

% 读取实部和虚部数据
y_ifft512_re = [];
y_ifft512_im = [];

% 读取实部数据
fid_re = fopen(re_file, 'r');
tline = fgets(fid_re);
while ischar(tline)
    y_ifft512_re = [y_ifft512_re; typecast(uint16(hex2dec(tline(1:4))), 'int16')];
    tline = fgets(fid_re);
end
fclose(fid_re);

% 读取虚部数据
fid_im = fopen(im_file, 'r');
tline = fgets(fid_im);
while ischar(tline)
    y_ifft512_im = [y_ifft512_im; typecast(uint16(hex2dec(tline(1:4))), 'int16')];
    tline = fgets(fid_im);
end
fclose(fid_im);

% 合并实部和虚部为复数数组
y_ifft = complex(double(y_ifft512_re) / (2^7), double(y_ifft512_im) / (2^7));

subplot(2, 3, 6); 
plot(abs(y_ifft), 'LineWidth', 1.5);
title('定点化RTL 512点IFFT结果（幅度）');
xlabel('频率索引');
ylabel('幅度');
grid on;
