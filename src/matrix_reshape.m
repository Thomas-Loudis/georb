function [data_matrix_2] = matrix_reshape(data_matrix, No_col_time, matrix_step)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  write_georb_format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Write coputations to output file according to the GEORB data format 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - data_matrix   : Matrix Data (Partials over one arc)
% - No_col_time   : Number of collums for time argument(s); 1 or 2 collumns e.g. MJD or MJD, Sec since 00 h
% - matrix_step   : Number of rows of the matrix per epoch (e.g. 6)
%
% Output arguments:
% - data_matrix_2 : Matrix Data reshaped (row per epoch approach)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                           14 February 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Reshape Partials matrix-based format to row per epoch format
[d1, d2] = size(data_matrix);
i_data2 = 0;
Nparam_col = d2 - No_col_time;
data_matrix_2 = zeros(d1/matrix_step, No_col_time + Nparam_col * matrix_step);
for i_data = 1 : matrix_step : d1
    i_data2 = i_data2 + 1;
    data_matrix_2(i_data2, 1:No_col_time) = data_matrix(i_data , 1 : No_col_time);    
    for i_step = 1 : matrix_step 
        col_1 = No_col_time + (i_step-1) * Nparam_col + 1 ;
        col_2 = No_col_time + (i_step) * Nparam_col;
        data_matrix_2(i_data2 , col_1 : col_2) = data_matrix(i_data + i_step-1 , No_col_time + 1 : d2);
    end
end
