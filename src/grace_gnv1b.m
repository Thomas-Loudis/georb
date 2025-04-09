function [gnv1b,COVmatrix] = grace_gnv1b(filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grace_kbr:  GRACE GNV1B orbit data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Reading GRACE satellite orbit from the ISDC data files GPS Navigation
%  Level 1B data (GNV1B).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - filename:   GNV1B data file's name
%
% Output arguments:
% - kbr1b:      KBR data array
%   gnv1b = [MJDgps tgps GNV1B_position GNV1B_veloctiy ]
%   MJDgps:     MJD in GPS time including fraction of the day
%   tgps:       Seconds since 0h in GPS time scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   April 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 10/08/2012   Covariance matrix has been formed and added to the output
%              arguments
% 11/05/2021, Dr. Thomas Papanikolaou
%             Modified to read both GRACE and GRACE-FO mission data format 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Covariance matrix formulation: Reading or not
COV_flag = 0;

% ISDC GRACE data GPS Time start epoch: 12:00 01-Jan-2000
[jd2000,mjd2000] = mjd3(12*3600,1,1,2000);

fid = fopen(filename);
i = 0;
j = 1;
Nepochs = 0;
while (~feof(fid))
    line = fgetl(fid);
    if i == 1
        Nepochs = Nepochs + 1;       
    end    
    % GRACE satellites (1 and 2 or A and B)
    str_endheader = sscanf(line,'%13c %*');
    endofheader_grace = 'END OF HEADER';
    test = strcmp(str_endheader,endofheader_grace);
    if test == 1
        i = 1;
       position_end_of_header = ftell(fid);             
    end    
    % Releases 2, 3
    str_endheader = sscanf(line,'%s%c%s%c%s%c%s%*');
    endofheader_grace = '# END OF HEADER';
    test = strcmp(str_endheader,endofheader_grace);
    if test == 1
        i = 1;
       position_end_of_header = ftell(fid);             
    end        
    % GRACE-FO satellites (3 and 4 or C and D)
    str_endheader = sscanf(line,'%20c %*');
    endofheader_gracefo = '# End of YAML header';
    test = strcmp(str_endheader,endofheader_gracefo);
    if test == 1
        i = 1;
       position_end_of_header = ftell(fid);             
    end
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of arrays
gnv1b = zeros(Nepochs,8);
if COV_flag > 0
COVmatrix_cv = zeros(3 * Nepochs, 3 * Nepochs); 
COVtime = zeros(3 * Nepochs, 1);
COVmatrix = zeros(3 * Nepochs, 3 * Nepochs + 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename);
fseek(fid,position_end_of_header,'bof');
i = 0;
j = 1;
while (~feof(fid))
    line = fgetl(fid);
    [data_ith_vec] = sscanf(line,'%d %*s %*s %e %e %e %e %e %e %e %e %e %e %e %e %*');
    tgps2000 = fix(data_ith_vec(1,1));
    mjd = mjd2000 + tgps2000 / (24*3600);
    tgps2000_o = (fix(mjd) - mjd2000) * (24*3600);
    % Seconds sice 0h
    tgps = tgps2000 - tgps2000_o;
    %[tgps,D,M,Y] = MJD_inv(mjd);
    % gnv1bdata = [str2num(line(15:end))];
    % gnv1b(j,:) = [mjd tgps gnv1bdata(1,1:3) gnv1bdata(1,7:9)];
    gnv1bdata_pos = data_ith_vec(2:4,1)';
    gnv1bdata_vel = data_ith_vec(8:10,1)';
    gnv1b(j,:) = [mjd tgps gnv1bdata_pos gnv1bdata_vel];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Covariace matrix of positions
if COV_flag > 0
        COVmatrix_cv(j*3-2 : j*3 , j*3-2 : j*3) = [
            gnv1bdata(1,4)^2            0                   0
                0               gnv1bdata(1,5)^2            0
                0                       0           gnv1bdata(1,6)^2];
        COVtime(j*3-2 : j*3 , 1) = [mjd mjd mjd]';    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j = j + 1;
end
fclose(fid);

if COV_flag > 0
    COVmatrix = [COVtime COVmatrix_cv];   
elseif COV_flag == 0
    COVmatrix = 0;
end
