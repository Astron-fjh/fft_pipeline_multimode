%% Description:
%   ������ʵ���˿��ٸ���Ҷ�任��FFT���еĵ������㣨Butterfly Operation����
%   �ú����������������źţ�xm1��xm2������ת���ӣ�w����
%   ͨ�����ι�ʽ���㲢��������µĸ����źţ�xn1��xn2����
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
    
    xn1_re = floor(xn1_re / (2^13));    % ������ת���������� 2^13����������� 2^13 ���л�ԭ
    xn1_im = floor(xn1_im / (2^13));
    
    xn2_re = floor(xn2_re / (2^13));
    xn2_im = floor(xn2_im / (2^13));
    
    xn1_re = int16(xn1_re);  % �ض�Ϊ 16 λ�з�������
    xn1_im = int16(xn1_im);

    xn2_re = int16(xn2_re);
    xn2_im = int16(xn2_im);
end