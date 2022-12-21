function [GMconstant,AU,EMRAT,deformat,deperiod] = dexxxeph_readhd(HDfilename,GMearth)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dexxxeph_readhd :  JPL's DExxx header file reading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read header file header.xxx of JPL's DE series planetary and lunar
%  ephemerides.
%
% Input arguments
% - HDfilename : Header file name of the DExxx ephemeris
% - GMearth    : Earth's GM gravity constant
%                If argument is missing the following value is used instead
%                GMearth = 3.986004415*10^14
%
% Output arguments:
% - GMconstant : DExxx GM gravity constants in m^3/sec^2
%   >> GMconstant = [GM1 GM2 GMB GM4 GM5 GM6 GM7 GM8 GM9 GMmoon GMS]'  
% - AU         : AU constant (Km/au)
% - EMRAT      : GMearth/GMmoon ratio
% - deformat   : DExxx data record format
% - deperiod   : DExxx data record period (e.g. 32 days)
%
%  GM gravity constants are converted from au^3/day^2 to m^3/sec^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin == 1
    % units : m^3/sec^2
    GMearth = 3.986004415 * 10^14;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read DExxx Header file: header.xxx
fid = fopen(HDfilename);
i = 0;
j = 0;
line_no = 0;
linei_group1030 = 1000;
linei_group1041 = 1000;
linei_group1050 = 1000;
while (~feof(fid))
    line = fgetl(fid);
    i = i + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NCOEFF
    if i == 1
        line1num = sscanf(line,'%s %d');
        [sz1 sz2] = size(line1num);
        ncoeff = line1num(sz1,sz2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    linei_group1030 = linei_group1030 + 1;
    if length(line) > 11
        if line(1:12) == 'GROUP   1030'
            linei_group1030 = 1;
        end
    end
    if linei_group1030 == 3
        group1030 = str2num(line);
        deperiod = group1030(1,3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    linei_group1041 = linei_group1041 + 1;
    if length(line) > 11
        if line(1:12) == 'GROUP   1041'
            linei_group1041 = 1;
        end
    end
    if linei_group1041 == 6
        group1041_line = str2num(line);
        AU = group1041_line(1,1);
        EMRAT = group1041_line(1,2);
        GM1 = group1041_line(1,3);
    end
    if linei_group1041 == 7
        group1041_line = str2num(line);
        GM2 = group1041_line(1,1);
        GMB = group1041_line(1,2);
        GM4 = group1041_line(1,3);
    end
    if linei_group1041 == 8
        group1041_line = str2num(line);
        GM5 = group1041_line(1,1);
        GM6 = group1041_line(1,2);
        GM7 = group1041_line(1,3);
    end
    if linei_group1041 == 9
        group1041_line = str2num(line);
        GM8 = group1041_line(1,1);
        GM9 = group1041_line(1,2);
        GMS = group1041_line(1,3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    linei_group1050 = linei_group1050 + 1;
    if length(line) > 11
        if line(1:12) == 'GROUP   1050'
            line_no = i;
            linei_group1050 = 1;
        end
    end
    if linei_group1050 == 3
        deformat(1,:) = str2num(line);
    end
    if linei_group1050 == 4
        deformat(2,:) = str2num(line);
    end
    if linei_group1050 == 5
        deformat(3,:) = str2num(line);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
fclose(fid);
clear i j fid sz1 sz2 line_no
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing the number of Chebychev coefficients for each of the 13 items
[sz1 sz2] = size(deformat);
for j = 1 : sz2
    if j == 12
        deformat(4,j) = 2 * deformat(2,j) * deformat(3,j);
    else
        deformat(4,j) = 3 * deformat(2,j) * deformat(3,j);
    end
end
clear sz1 sz2 j 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DExxx constatnts: AU EMRAT GM1 GM2 GMB GM4 GM5 GM6 GM7 GM8 GM9 GMmoon GMS
% GMmoon in m^3/sec^2
GMmoon = GMearth / EMRAT;
GMconstant_au = [GM1 GM2 GMB GM4 GM5 GM6 GM7 GM8 GM9 GMmoon GMS]';
GMconstant = zeros(11,1);
for i = 1 : 11
    if i == 10
        GMconstant(i,1) = GMconstant_au(i,1);
    else
        % GM units conversion from au^3/day^2 to to m^3/sec^2
        GMconstant(i,1) = GMconstant_au(i,1) * (10^9 * AU^3) / (24 * 3600)^2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
