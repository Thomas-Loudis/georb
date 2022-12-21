function [quat] = grace_sca1b(filename,tstop)


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


% Reading or not the whole day data
if tstop == 0
    tstop = 24 * 3600;
end

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
    end
    % GRACE-FO satellites (3 and 4 or C and D)
    str_endheader = sscanf(line,'%20c %*');
    endofheader_gracefo = '# End of YAML header';
    test = strcmp(str_endheader,endofheader_gracefo);
    if test == 1
        i = 1;
    end
end
fclose(fid);
clear i j fid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of arrays
quat = zeros(Nepochs,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename);
% Read ISDC file
i = 0;
j = 1;
while (~feof(fid))
    line = fgetl(fid);
    if i == 1
        tgps2000 = str2num(line(1:9));
        mjd = mjd2000 + tgps2000 / (24*3600);
        tgps2000_o = (fix(mjd) - mjd2000) * (24*3600);
        % Seconds sice 0h
        tgps = tgps2000 - tgps2000_o;
        if tgps > tstop
            break
        end
        sca = [str2num(line(15:end))];
        quat(j,:) = [mjd tgps sca(1,1:4)];
        j = j + 1;       
        clear tgps mjd sca
    end
    % GRACE satellites (1 and 2 or A and B)
    str_endheader = sscanf(line,'%13c %*');
    endofheader_grace = 'END OF HEADER';
    test = strcmp(str_endheader,endofheader_grace);
    if test == 1
        i = 1;
    end    
    % GRACE-FO satellites (3 and 4 or C and D)
    str_endheader = sscanf(line,'%20c %*');
    endofheader_gracefo = '# End of YAML header';
    test = strcmp(str_endheader,endofheader_gracefo);
    if test == 1
        i = 1;
    end    
end
fclose(fid);
clear i j fid
