%% Description:
%   用于生成定点化的旋转因子，量化位宽为16bits
%   Used to generate fixed-point rotation factors with a quantization bit-width of 16 bits.
%
%% Author(s):
%   Astron-fjh

clear all;
clc;
N = 512;
for m = 0 : N/2-1
    W(m+1)=complex(cos(2*pi*m/N), -sin(2*pi*m/N));
end

re = real(W);
re = floor(re'*(2^13));
im = imag(W)
im = floor(im'*(2^13));

save('wn_re_512_fixed16.mat', 're');  % 保存实部
save('wn_im_512_fixed16.mat', 'im');  % 保存虚部

fidr = fopen('wn_re_512_fixed16.coe', 'wt'); 
fidi = fopen('wn_im_512_fixed16.coe', 'wt');

%- standard format
fprintf(fidr, 'MEMORY_INITIALIZATION_RADIX = 10;\n');                     
fprintf(fidr, 'MEMORY_INITIALIZATION_VECTOR =\n');
fprintf(fidi, 'MEMORY_INITIALIZATION_RADIX = 10;\n');                     
fprintf(fidi, 'MEMORY_INITIALIZATION_VECTOR =\n');

%- write data in coe file
for i = 1 : 1 : N/2
    fprintf(fidr, '%d,\n', re(i));  
    fprintf(fidi, '%d,\n', im(i));  
end

fclose(fidr);
fclose(fidi);