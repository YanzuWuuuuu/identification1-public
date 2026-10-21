function [output_matrix] = diag_real_imag_divide(input_matrix)
%將一矩陣的實部虛部分開，將其对角排列， []-> diag([real,imag])
%   輸入m*n輸出2m*2n
real_of_input = real(input_matrix);
imag_of_input = imag(input_matrix);
row_length = length(real_of_input(:,1));
colunm_length = length(real_of_input(1,:));
output_matrix = zeros(2*row_length,2*colunm_length);
output_matrix(1:row_length,1:colunm_length) = real_of_input;
output_matrix(row_length+1:end,colunm_length+1:end) = imag_of_input;

end

