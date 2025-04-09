function [quat] = grace_sca1b(filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grace_sca1b:  GRACE star camera assembly data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Reading GRACE star camera assembly data from the ISDC data files SCA1b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - filename:      SCA1b data file's name
%
% Output arguments:
% - sca:          Star Camera data: Quaternions
%   sca = [MJDgps tgps quatangle quanticoeff quatjcoeff quatkcoeff]
%   MJDgps:        MJD in GPS time including fraction of the day
%   tgps:          Seconds since 0h in GPS time scale
%   quatangle:     Cos mu/2 element of quartenion
%   quaticoeff:    I element of quaternion rotation axis
%   quatjcoeff:    J element of quaternion rotation axis
%   quatkcoeff:    K element of quaternion rotation axis
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
    str_line = sscanf(line,'%s');
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
quat = zeros(Nepochs,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename);
fseek(fid,position_end_of_header,'bof');
% Read ISDC file
i = 0;
j = 1;
while (~feof(fid))
    line = fgetl(fid);
        [data_ith_vec] = sscanf(line,'%d %*s %d %e %e %e %e %e %*s %*');        
        % tgps2000 = str2num(line(1:9));
        tgps2000 = fix(data_ith_vec(1,1));
        mjd = mjd2000 + tgps2000 / (24*3600);
        tgps2000_o = (fix(mjd) - mjd2000) * (24*3600);
        % Seconds sice 0h
        tgps = tgps2000 - tgps2000_o;
        % sca = [str2num(line(15:end))];
        sca = data_ith_vec(3:7,1)';
        quat(j,:) = [mjd tgps sca(1,1:4)];
        j = j + 1;       
end
fclose(fid);
