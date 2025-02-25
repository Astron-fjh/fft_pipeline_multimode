%% Description:
%   ������ʵ���˿��ٸ���Ҷ�任��FFT���еĻ�2�������㣨Radix-2 Butterfly Operation����
%   �ú�������һ��ʱ���ź� `xn`����ͨ���ݹ�ĵ�������������Ƶ���ʾ��FFT����
%   �����źŵĳ��Ȼᱻ����Ϊ��ӽ���2���ݴΣ��Ӷ�����FFT�㷨��Ҫ��
%
%   This function implements the Radix-2 Butterfly Operation in the Fast Fourier Transform (FFT). 
%   The function takes a time-domain signal xn and computes its frequency-domain representation (FFT) using recursive butterfly operations. 
%   The length of the input signal will be padded to the nearest power of 2 to meet the requirements of the FFT algorithm. 
%   This algorithm implements the classic Cooley-Tukey FFT method, utilizing bit reversal and butterfly operations to improve computational efficiency.
%   
%% Author(s):
%   Astron-fjh

function [Xk] = fft_radix2(xn, N)
    % ���������ź� xn �ĳ��� N ����Ҫ����С2���ݴ�
    M = nextpow2(N);  
    N_padded = 2^M;

    if length(xn) < N_padded
        xn = [xn, zeros(1, N_padded - length(xn))];  % ����
    end

    n = 1:N_padded;
    x = xn(bitrevorder(n-1) + 1);   % λ��ת

    % ���ļ��м���512�㶨�㻯��ת����
    load('wn_re_512_fixed16.mat', 're');
    load('wn_im_512_fixed16.mat', 'im');
    
    W = complex(re, im);

    % FFT����
    for m = 1:M            % mΪ����
        for i = 1:2^(M-m)        % iΪ����
            for j = 1:2^(m-1)         % jΪ���������ڵı��
                k = (i-1) * 2^m + j;
                n_idx = 2^(M-m)*(j-1);  % ���ڲ�ͬ�� N����������

                % ���� N ��������ת���ӵ�����
                if N == 256
                    W_idx = n_idx * 2;    % ÿ������ת����ȡһ��
                elseif N == 128
                    W_idx = n_idx * 4;    % ÿ�ĸ���ת����ȡһ��
                elseif N == 64
                    W_idx = n_idx * 8;    % ÿ�˸���ת����ȡһ��
                else
                    W_idx = n_idx;        % ���� 512�㣬ֱ��ʹ��
                end
                
                [T1_re, T1_im, T2_re, T2_im] = butterfly(real(x(k)), imag(x(k)), ...
                                                         real(x(k+2^(m-1))), imag(x(k+2^(m-1))), ...
                                                         real(W(W_idx+1)), imag(W(W_idx+1)));
                
                if m == 1  
                    % ��ӡ��ǰ k �� k+2^(m-1) ��ֵ�Լ���ת���� W(W_idx+1)
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
        
        % ÿ������һ������� x ��ֵ
        % disp(['After level ', num2str(m), ' computation:']);
        % disp(x);
    end

    Xk = x;
end
