function [orbit_data_matrix, header_data_matrix] = read_georb_data(data_filename)


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
% Last modified
% 8 March 2025 Thomas Loudis
%              Optimizing code for reading big data, minimizing CPU time by
%              using only fscanf and removing the use of fgetl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(data_filename);
data_rec_start = 0;
N_data_rec = 0;
i_header_rec = 0;
while (~feof(fid))
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');
    if data_rec_start == 1
        % Next data record
        N_data_rec = N_data_rec + 1;
        [data_rec_part, N_data_entry]  = sscanf(line_ith,'%s');
        break
    end    
    test = strcmp(str1,'end_of_header');
    if test == 1
       data_rec_start = 1;
       position_end_of_header = ftell(fid);     
    end
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix preallocation
Nelements = N_data_entry;
% orbit_data_matrix = zeros(N_data_rec, Nelements);

% header_data_matrix(1 : i_header_rec).header_data = ' ';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(data_filename);
data_rec_start = 0;
N_data_rec = 0;
i_header_rec = 0;
while (~feof(fid))
% while (data_rec_start < 1)
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');
    % Read Header
    if data_rec_start == 0
        i_header_rec = i_header_rec + 1;
        % header_data_matrix.header_data(i_header_rec,:) = line_ith;
        header_data_matrix(i_header_rec).header_data = line_ith;
        % header_line_read_str1 = header_data_matrix(i_header_rec).header_data
        % header_line_str1 = header_data_matrix(i_header_rec).header_data
    end
    % End of Header test 
    test = strcmp(str1,'end_of_header');
    if test == 1
       data_rec_start = 1;
       break
    end    
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(data_filename);
fseek(fid,position_end_of_header,'bof');
size_data_matrix = [Nelements Inf];
% '%d %d %e'
% formatSpec = '%e \n'
formatSpec = '%e';
data_matrix = fscanf(fid,formatSpec,size_data_matrix);
fclose(fid);
orbit_data_matrix = data_matrix';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
