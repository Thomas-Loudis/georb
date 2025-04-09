function [acc] = grace_acc1b(filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grace_acc1b:  GRACE accelerometry data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Reading GRACE accelerometry data from the ISDC data files ACC1b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - filename:      accelerometry data file's name
%
% Output arguments:
% - acc:           Accelerometry array
%   acc = [MJDgps tgps lin_accl_x lin_accl_y lin_accl_z]
%   MJDgps:        MJD in GPS time including fraction of the day
%   tgps:          Seconds since 0h in GPS time scale
%   lin_accl_i:    Linear acceleration along Science Reference Frame axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   April 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 11/05/2021, Dr. Thomas Papanikolaou
%             Modified to read both GRACE and GRACE-FO mission data format 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
    % Release 1
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
acc = zeros(Nepochs,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename);
fseek(fid,position_end_of_header,'bof');
% Read ISDC file
i = 0;
j = 1;
while (~feof(fid))
    line = fgetl(fid);
    [data_ith_vec] = sscanf(line,'%d %*s %e %e %e %*');
    % tgps2000 = str2num(line(1:9));
    tgps2000 = fix(data_ith_vec(1,1));
    mjd = mjd2000 + tgps2000 / (24*3600);
    tgps2000_o = (fix(mjd) - mjd2000) * (24*3600);
    % Seconds sice 0h
    tgps = tgps2000 - tgps2000_o;
    %[tgps,D,M,Y] = MJD_inv(mjd);
    % accl = [str2num(line(13:end))];
    accl = data_ith_vec(2:4,1)';
    acc(j,:) = [mjd tgps accl(1,1:3)];
    j = j + 1;
end
fclose(fid);
