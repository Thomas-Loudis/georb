function [orbit_data_matrix] = read_georb_data_v2(data_filename) 

%-------------------------------------------------------------------------- 
% Copyright (C) 2007-present  Thomas (Loudis) Papanikolaou  
%  
% GEORB :: 'Gravity and Precise Orbit Determination system' 
%  
% GEORB is a software for Precise Orbit Determination (POD), gravity field  
% recovery and mission design 
%  
% GEORB is licensed under the GNU Affero General Public License v3.0 while  
% for information regarding dual-license options visit the official website 
% www.georb.gr  
%-------------------------------------------------------------------------- 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Function:  read_georb_data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  Read orbit data products files in georb format  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - data_filename       : Data file name in georb data format  
% 
% Output arguments: 
% - orbit_data_matrix   : Data matrix in format: [MJD Sec_00 orbit(i,:)] 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                           27 November 2022 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
fid = fopen(data_filename); 
data_rec_start = 0; 
N_data_rec = 0; 
while (~feof(fid)) 
    line_ith = fgetl(fid); 
    str1 = sscanf(line_ith,'%s %*');     
    if data_rec_start > 0 
        % Next data record 
        N_data_rec = N_data_rec + 1; 
        [data_rec_part, N_data_entry]  = sscanf(line_ith,'%s'); 
    end     
    test = strcmp(str1,'end_of_header'); 
    if test == 1 
       data_rec_start = 1; 
    end 
end 
fclose(fid); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Matrix preallocation 
Nelements = N_data_entry; 
orbit_data_matrix = zeros(N_data_rec, Nelements); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
fid = fopen(data_filename); 
data_rec_start = 0; 
N_data_rec = 0; 
while (~feof(fid)) 
    line_ith = fgetl(fid); 
    str1 = sscanf(line_ith,'%s %*'); 
    % Read data records 
    if data_rec_start > 0 
        % Next data record 
        N_data_rec = N_data_rec + 1; 
        % MJD number 
        % mjd = sscanf(line_ith,'%d %*'); 
        mjd = sscanf(line_ith,'%e %*'); 
        % Seconds since start of the day (00h) 
        % sec_00 = sscanf(line_ith,'%*d %f %*'); 
        sec_00 = sscanf(line_ith,'%*e %e %*'); 
        % Store Epoch's arguments to the orbit data matrix 
        orbit_data_matrix(N_data_rec,1) = mjd; 
        orbit_data_matrix(N_data_rec,2) = sec_00; 
        % Number of data elements 
        [data_rec_part, N_data_entry]  = sscanf(line_ith,'%s'); 
        % Get read index position after Epoch arguments 
        [data_rec_fields12,n,errmsg, read_index] = sscanf(line_ith,'%s%c',4); 
        Nchar_data = length(line_ith); 
        data_rec_part = line_ith(read_index : Nchar_data); 
        % Read orbit data argument but epoch's ones 
        orbit_data = sscanf(data_rec_part,'%e%*c', (N_data_entry-2)*2); 
        orbit_data_matrix(N_data_rec, 3:N_data_entry) = orbit_data'; 
    end 
    % End of Header test  
    test = strcmp(str1,'end_of_header'); 
    if test == 1 
       data_rec_start = 1; 
    end     
end 
fclose(fid); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
