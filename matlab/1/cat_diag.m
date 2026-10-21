function [output_matrix] = cat_diag(input_matrix1,input_matrix2)
[row_matrix1,col_matrix1] = size(input_matrix1);
[row_matrix2,col_matrix2] = size(input_matrix2);
output_matrix = zeros(row_matrix1+row_matrix2,col_matrix1+col_matrix2);
output_matrix(1:row_matrix1,1:col_matrix1) = input_matrix1;
output_matrix(row_matrix1+1:row_matrix1+row_matrix2,col_matrix1+1:col_matrix1+col_matrix2) = input_matrix2;
end

