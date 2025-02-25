%% Description:
%   本函数实现了快速傅里叶变换（FFT）中的蝶形运算（Butterfly Operation）。
%   该函数接受两个复数信号（xm1和xm2）及旋转因子（w），
%   通过蝶形公式计算并输出两个新的复数信号（xn1和xn2）。
%
%   This function implements the Butterfly Operation in the Fast Fourier Transform (FFT). 
%   The function takes two complex signals (xm1 and xm2) and a rotation factor (w), 
%   and computes and outputs two new complex signals (xn1 and xn2) using the butterfly formula.
%
%% Author(s):
%   Astron-fjh

function [xn1_re, xn1_im, xn2_re, xn2_im] = butterfly(xm1_re, xm1_im, xm2_re, xm2_im, w_re, w_im)
    xn1_re = xm1_re * (2^13) + xm2_re*w_re - (xm2_im*w_im);
    xn1_im = xm1_im * (2^13) + xm2_re*w_im + xm2_im*w_re;

    xn2_re = xm1_re * (2^13) - (xm2_re*w_re) + xm2_im*w_im;
    xn2_im = xm1_im * (2^13) - (xm2_re*w_im) - (xm2_im*w_re);
    
    xn1_re = floor(xn1_re / (2^13));    % 由于旋转因子扩大了 2^13，在这里除以 2^13 进行还原
    xn1_im = floor(xn1_im / (2^13));
    
    xn2_re = floor(xn2_re / (2^13));
    xn2_im = floor(xn2_im / (2^13));
    
    xn1_re = int16(xn1_re);  % 截断为 16 位有符号整数
    xn1_im = int16(xn1_im);

    xn2_re = int16(xn2_re);
    xn2_im = int16(xn2_im);
end