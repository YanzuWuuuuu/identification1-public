function [out] = real_ima_divide(C)
%REAL_IMA_DIVIDE 分开复数矩阵的实虚部，使其上下排列  [C] -> [R;I]
% input : C:complex matrix 
% output: real matrix [Real;Imag]
real_part = real(C);
imag_part = imag(C);
out = [real_part;imag_part];
end

