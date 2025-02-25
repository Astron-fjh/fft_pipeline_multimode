%% Description:
%   本函数实现了快速傅里叶变换（FFT）中的基2蝶形运算（Radix-2 Butterfly Operation）。
%   该函数接受一个时域信号 `xn`，并通过递归的蝶形运算计算出其频域表示（FFT）。
%   输入信号的长度会被补齐为最接近的2的幂次，从而适配FFT算法的要求。
%
%   This function implements the Radix-2 Butterfly Operation in the Fast Fourier Transform (FFT). 
%   The function takes a time-domain signal xn and computes its frequency-domain representation (FFT) using recursive butterfly operations. 
%   The length of the input signal will be padded to the nearest power of 2 to meet the requirements of the FFT algorithm. 
%   This algorithm implements the classic Cooley-Tukey FFT method, utilizing bit reversal and butterfly operations to improve computational efficiency.
%   
%% Author(s):
%   Astron-fjh

function [Xk] = fft_radix2(xn, N)
    % 计算输入信号 xn 的长度 N 所需要的最小2的幂次
    M = nextpow2(N);  
    N_padded = 2^M;

    if length(xn) < N_padded
        xn = [xn, zeros(1, N_padded - length(xn))];  % 补零
    end

    n = 1:N_padded;
    x = xn(bitrevorder(n-1) + 1);   % 位反转

    % 从文件中加载512点定点化旋转因子
    load('wn_re_512_fixed16.mat', 're');
    load('wn_im_512_fixed16.mat', 'im');
    
    W = complex(re, im);

    % FFT计算
    for m = 1:M            % m为级数
        for i = 1:2^(M-m)        % i为组数
            for j = 1:2^(m-1)         % j为蝶形在组内的编号
                k = (i-1) * 2^m + j;
                n_idx = 2^(M-m)*(j-1);  % 对于不同的 N，调整索引

                % 根据 N 来调整旋转因子的索引
                if N == 256
                    W_idx = n_idx * 2;    % 每两个旋转因子取一个
                elseif N == 128
                    W_idx = n_idx * 4;    % 每四个旋转因子取一个
                elseif N == 64
                    W_idx = n_idx * 8;    % 每八个旋转因子取一个
                else
                    W_idx = n_idx;        % 对于 512点，直接使用
                end
                
                [T1_re, T1_im, T2_re, T2_im] = butterfly(real(x(k)), imag(x(k)), ...
                                                         real(x(k+2^(m-1))), imag(x(k+2^(m-1))), ...
                                                         real(W(W_idx+1)), imag(W(W_idx+1)));
                
                if m == 1  
                    % 打印当前 k 和 k+2^(m-1) 的值以及旋转因子 W(W_idx+1)
                    disp(['x(k): ', num2str(real(x(k))), ' + ', num2str(imag(x(k))), 'i']);
                    disp(['x(k+2^(m-1)): ', num2str(real(x(k+2^(m-1)))), ' + ', num2str(imag(x(k+2^(m-1)))) , 'i']);
                    disp(['W(W_idx+1): ', num2str(real(W(W_idx+1))), ' + ', num2str(imag(W(W_idx+1))), 'i']);
                end

                x(k+2^(m-1)) = complex(T2_re, T2_im);
                x(k) = complex(T1_re, T1_im);
                
                if m == 1
                    disp(['After level ', num2str(m), ' computation for k = ', num2str(k)]);
                    disp(['x(k): ', num2str(real(x(k))), ' + ', num2str(imag(x(k))), 'i']);
                    disp(['x(k+2^(m-1)): ', num2str(real(x(k+2^(m-1)))), ' + ', num2str(imag(x(k+2^(m-1)))) , 'i']);
                end
            end
        end
        
        % 每计算完一级，输出 x 的值
        % disp(['After level ', num2str(m), ' computation:']);
        % disp(x);
    end

    Xk = x;
end
