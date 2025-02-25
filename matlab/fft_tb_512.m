%% Description:
%   ������ʵ����512��FFT��IFFT�Ķ��㻯���������׼FFT��IFFT���жԱȡ�
%   ���ȣ�����һ��ʱ���źŲ�����FFT��֮���ź�ת��Ϊ16λ�����ʾ������Ϊ�ı��ļ�������RTL�����롣
%   ���ʹ���Զ���Ķ��㻯FFT�㷨���м��㣬����������ͬʱ�������׼IFFT�����ж��㻯����
%   ��󣬴�����ӻ���׼FFT�����㻯FFT����׼IFFT�Ͷ��㻯IFFT�Ľ��������RTL������жԱȡ�
%
%   This code implements the fixed-point processing of 512-point FFT and IFFT, and compares the results with standard FFT and IFFT.
%   First, a time-domain signal is generated and processed with FFT. The signal is then converted to a 16-bit fixed-point representation and saved as a text file to be used as input for RTL.
%   Next, the custom fixed-point FFT algorithm is used for computation, and the results are saved. The standard IFFT is also computed and processed with fixed-point representation.
%   Finally, the code visualizes the results of the standard FFT, fixed-point FFT, standard IFFT, and fixed-point IFFT, and compares them with the RTL output.

%% Author(s):
%   Astron-fjh

% ��ʼ��
clear all
clc
N = 512;    % 512��
t = 1:1:N;
X_re = cos(1/3*pi*t);
X_im = zeros(1, N);
X = complex(X_re, X_im);
X_fft = (fft(X));

% ʱ���ź�ת��Ϊ�����Ʋ��벢����
X_re_16bit = floor(X_re * (2^7));   % ���㻯Ϊ16λ�������Ŵ��� 2^7
X_im_16bit = floor(X_im * (2^7));

% ��֤ת��Ϊ�з���16λ����
X_re_int16 = int16(X_re_16bit);     % ת��Ϊ�з���16λ����
X_im_int16 = int16(X_im_16bit);

% ���ļ���д������
fd_re = fopen("./512��/x_re512.txt", 'wb');
fd_im = fopen("./512��/x_im512.txt", 'wb');
for i = 1:N
    % ��16λ�з�������תΪ���벢����
    fprintf(fd_re, '%04x\r\n', typecast(X_re_int16(i), 'uint16'));
    fprintf(fd_im, '%04x\r\n', typecast(X_im_int16(i), 'uint16'));
end
fclose(fd_re);
fclose(fd_im);

% Ƶ���źţ�FFT�����ת��Ϊ�����Ʋ��벢����
X_fft_re_16bit = floor(real(X_fft) * (2^7));  % �Ŵ���2^7
X_fft_im_16bit = floor(imag(X_fft) * (2^7));

% ת��Ϊ16λ�з�������
X_fft_re_int16 = int16(X_fft_re_16bit);
X_fft_im_int16 = int16(X_fft_im_16bit);

fd_re = fopen("./512��/x_fft512_re.txt", 'wb');
fd_im = fopen("./512��/x_fft512_im.txt", 'wb');
for i = 1:N
    fprintf(fd_re, '%04x\r\n', typecast(X_fft_re_int16(i), 'uint16'));
    fprintf(fd_im, '%04x\r\n', typecast(X_fft_im_int16(i), 'uint16'));
end
fclose(fd_re);
fclose(fd_im);

% ���㻯FFT����
X_16bit = complex(X_re_16bit, X_im_16bit);
X_sim_fft = fft_radix2(X_16bit, N);

% ת�������涨�㻯FFT������
X_sim_fft_re_16bit = floor(real(X_sim_fft));
X_sim_fft_im_16bit = floor(imag(X_sim_fft));

X_sim_fft_re_int16 = int16(X_sim_fft_re_16bit);
X_sim_fft_im_int16 = int16(X_sim_fft_im_16bit);

fd_re = fopen("./512��/x_sim_fft512_re.txt", 'wb');
fd_im = fopen("./512��/x_sim_fft512_im.txt", 'wb');
for i = 1:N
    fprintf(fd_re, '%04x\r\n', typecast(X_sim_fft_re_int16(i), 'uint16'));
    fprintf(fd_im, '%04x\r\n', typecast(X_sim_fft_im_int16(i), 'uint16'));
end
fclose(fd_re);
fclose(fd_im);

X_sim_fft_re = floor(real(X_sim_fft)) / (2^7);
X_sim_fft_im = floor(imag(X_sim_fft)) / (2^7);
X_sim_fft = complex(X_sim_fft_re, X_sim_fft_im);

% ���ӻ��Ƚ�
figure;
set(gcf, 'Position', [100, 100, 1200, 900]);
subplot(2, 3, 1); 
plot(abs(X_fft), 'LineWidth', 1.5);
title('��׼512��FFT��������ȣ�');
xlabel('Ƶ������');
ylabel('����');
grid on;

subplot(2, 3, 2); 
plot(abs(X_sim_fft), 'LineWidth', 1.5);
title('���㻯512��FFT��������ȣ�');
xlabel('Ƶ������');
ylabel('����');
grid on;

%%
% ��ȡRTL�����
re_file = './512��/y_fft512_re.txt';
im_file = './512��/y_fft512_im.txt';

% ��ȡʵ�����鲿����
y_fft512_re = [];
y_fft512_im = [];

% ��ȡʵ������
fid_re = fopen(re_file, 'r');
tline = fgets(fid_re);
while ischar(tline)
    y_fft512_re = [y_fft512_re; typecast(uint16(hex2dec(tline(1:4))), 'int16')];
    tline = fgets(fid_re);
end
fclose(fid_re);

% ��ȡ�鲿����
fid_im = fopen(im_file, 'r');
tline = fgets(fid_im);
while ischar(tline)
    y_fft512_im = [y_fft512_im; typecast(uint16(hex2dec(tline(1:4))), 'int16')];
    tline = fgets(fid_im);
end
fclose(fid_im);

% �ϲ�ʵ�����鲿Ϊ��������
y_fft = complex(double(y_fft512_re) / (2^7), double(y_fft512_im) / (2^7));

subplot(2, 3, 3); 
plot(abs(y_fft), 'LineWidth', 1.5);
title('���㻯RTL 512��FFT��������ȣ�');
xlabel('Ƶ������');
ylabel('����');
grid on;

%%
% ��׼IFFT����
X_ifft = ifft(X);

% ��׼IFFT�Ķ��㻯������FFT��ͬ��
X_ifft_re_16bit = floor(real(X_ifft) * (2^7));  % �Ŵ���2^7
X_ifft_im_16bit = floor(imag(X_ifft) * (2^7));  % �Ŵ���2^7

% ת��Ϊ16λ�з�������
X_ifft_re_int16 = int16(X_ifft_re_16bit);
X_ifft_im_int16 = int16(X_ifft_im_16bit);

% �����׼IFFT�Ķ��㻯���
fd_re = fopen("./512��/x_ifft512_re.txt", 'wb');
fd_im = fopen("./512��/x_ifft512_im.txt", 'wb');
for i = 1:N
    fprintf(fd_re, '%04x\r\n', typecast(X_ifft_re_int16(i), 'uint16'));
    fprintf(fd_im, '%04x\r\n', typecast(X_ifft_im_int16(i), 'uint16'));
end
fclose(fd_re);
fclose(fd_im);

% ���㻯IFFT����
% ����IFFT���Ƚ�FFT����������ȡ����ٶ�FFT���ȡ������ڵ���������������6λ
X_16bit_ifft = complex(X_re_16bit, X_im_16bit);
X_16bit_ifft = conj(X_16bit_ifft);  % ȡ����
X_sim_ifft = conj(fft_radix2(X_16bit_ifft, N)) / (2^9); 

% ���涨�㻯IFFT���
X_sim_ifft_re_16bit = floor(real(X_sim_ifft));
X_sim_ifft_im_16bit = floor(imag(X_sim_ifft));

X_sim_ifft_re_int16 = int16(X_sim_ifft_re_16bit);
X_sim_ifft_im_int16 = int16(X_sim_ifft_im_16bit);

fd_re = fopen("./512��/x_sim_ifft512_re.txt", 'wb');
fd_im = fopen("./512��/x_sim_ifft512_im.txt", 'wb');
for i = 1:N
    fprintf(fd_re, '%04x\r\n', typecast(X_sim_ifft_re_int16(i), 'uint16'));
    fprintf(fd_im, '%04x\r\n', typecast(X_sim_ifft_im_int16(i), 'uint16'));
end
fclose(fd_re);
fclose(fd_im);

X_sim_ifft = floor(X_sim_ifft) / (2^7);

% ���ӻ���׼IFFT�붨�㻯IFFT�Ա�
subplot(2, 3, 4); 
plot(abs(X_ifft), 'LineWidth', 1.5);
title('��׼512��IFFT��������ȣ�');
xlabel('Ƶ������');
ylabel('����');
grid on;

subplot(2, 3, 5); 
plot(abs(X_sim_ifft), 'LineWidth', 1.5);
title('���㻯512��IFFT��������ȣ�');
xlabel('Ƶ������');
ylabel('����');
grid on;

%%
% ��ȡRTL�����
re_file = './512��/y_ifft512_re.txt';
im_file = './512��/y_ifft512_im.txt';

% ��ȡʵ�����鲿����
y_ifft512_re = [];
y_ifft512_im = [];

% ��ȡʵ������
fid_re = fopen(re_file, 'r');
tline = fgets(fid_re);
while ischar(tline)
    y_ifft512_re = [y_ifft512_re; typecast(uint16(hex2dec(tline(1:4))), 'int16')];
    tline = fgets(fid_re);
end
fclose(fid_re);

% ��ȡ�鲿����
fid_im = fopen(im_file, 'r');
tline = fgets(fid_im);
while ischar(tline)
    y_ifft512_im = [y_ifft512_im; typecast(uint16(hex2dec(tline(1:4))), 'int16')];
    tline = fgets(fid_im);
end
fclose(fid_im);

% �ϲ�ʵ�����鲿Ϊ��������
y_ifft = complex(double(y_ifft512_re) / (2^7), double(y_ifft512_im) / (2^7));

subplot(2, 3, 6); 
plot(abs(y_ifft), 'LineWidth', 1.5);
title('���㻯RTL 512��IFFT��������ȣ�');
xlabel('Ƶ������');
ylabel('����');
grid on;
