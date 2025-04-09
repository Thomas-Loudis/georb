function [kbr1b,biasrange,rangerate,rangeaccl] = grace_kbr1b(filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grace_kbr:  GRACE K-band range (KBR) observations data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Reading GRACE KBR data from the ISDC data files KBR1b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - filename:   KBR1b data file's name
%
% Output arguments:
% - kbr1b:      KBR data array
%   kbr1b = [MJDgps tgps KBR1B(2:end)]
% - biasrange:  Corrected biased range data
%   biasrange = [MJDgps tgps biased_range]
% - rangerate:  Corrected range rate data
%   rangerate = [MJDgps tgps range_rate]
% - rangeaccl:  Corrected range acceleration data
%   rangeaccl = [MJDgps tgps range_acceleration]
%
%   MJDgps:     MJD in GPS time including fraction of the day
%   tgps:       Seconds since 0h in GPS time scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   April 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 12/05/2021, Dr. Thomas Papanikolaou
%             Modified to read both GRACE and GRACE-FO mission data format 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ISDC GRACE data GPS Time start epoch: 12:00 01-Jan-2000
[jd2000,mjd2000] = mjd3(12*3600,1,1,2000);

% Compute the number of epochs of the data
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
    %if i == 1
    if Nepochs == 1
       line = fgetl(fid);
       kbrdata = [str2num(line(1:end))];
       [N1, Nelements] = size(kbrdata);
    end
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of arrays
kbr1b = zeros(Nepochs,Nelements+1);
biasrange = zeros(Nepochs,3);
rangerate = zeros(Nepochs,3);
rangeaccl = zeros(Nepochs,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read and store data
fid = fopen(filename);
fseek(fid,position_end_of_header,'bof');
% Read ISDC file
i = 0;
j = 1;
while (~feof(fid))
    line = fgetl(fid);
%     if i == 1
        tgps2000 = str2num(line(1:9));
        mjd = mjd2000 + tgps2000 / (24*3600);
        tgps2000_o = (fix(mjd) - mjd2000) * (24*3600);
        % Seconds sice 0h
        tgps = tgps2000 - tgps2000_o;
        %[tgps,D,M,Y] = MJD_inv(mjd);
        kbrdata = [str2num(line(1:end))];
        kbr1b(j,:) = [mjd tgps kbrdata(2:end)];
        % Corrected Biased Range
        biased_range_corrected = kbrdata(1,2) + kbrdata(1,6) + kbrdata(1,9);
        biasrange(j,:) = [mjd tgps biased_range_corrected];
        % Corrected Range Rate
        range_rate_corrected = kbrdata(1,3) + kbrdata(1,7) + kbrdata(1,10);
        rangerate(j,:) = [mjd tgps range_rate_corrected];
        % Corrected Range Acceleration
        range_accl_corrected = kbrdata(1,4) + kbrdata(1,8) + kbrdata(1,11);
        rangeaccl(j,:) = [mjd tgps range_accl_corrected];       
        j = j + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Phase-break epochs
% 1 June 2021: Thomas Papanikolaou, Revision for GRACE Follow-on data
        qualflg = sscanf(line,'%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %s %*');
        if qualflg(1,1) == '1'
            phasebreak_epoch = [mjd tgps];
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end
fclose(fid);
