function [decbv,datarecord,deformat] = dexxxeph_read(DEfilename,HDfilename,tjd)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dexxxeph_read :  JPL's DExxx Planetary and Lunar ephemerides reading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read JPL's DE series planetary and lunar ephemerides.
%  Read the required data blocks according to the time input arguments and
%  store the Chebyshev coefficients of the respective data record for all
%  the 13 items (solar system  bodies, nutations and librations).
%
% Input arguments
% - DEfilename : DExxx file name
% - HDfilename : Header file name of the DExxx ephemeris
% - tjd        : Julian Date number of the required epoch
%
% Output arguments:
% - decbv      : DExxx Chebyshev coefficients for the required data record
% - datarecord : DExxx data record's number, start, end JD time
%   >>           datarecord = [ datarecord_no JDstart JDend]'
% - deformat   : DExxx data record format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read: Header file
fid = fopen(HDfilename);
i = 0;
j = 0;
line_no = 0;
while (~feof(fid))
    line = fgetl(fid);
    i = i + 1;
    if i == 1
        line1num = sscanf(line,'%s %d');
        [sz1 sz2] = size(line1num);
        ncoeff = line1num(sz1,sz2);
    end
    if length(line) > 11
        if line(1:12) == 'GROUP   1050'
            line_no = i;
        end
    end
    if line_no > 0
        if i >= (line_no+2) & i <= (line_no+4)
            j = j + 1;
            deformat(j,:) = str2num(line);
        end
    end
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
% Read: DExxx file
fid = fopen(DEfilename);
% Index of line
linei = 0;
% Index of data block
dblinei = -1000000;
while (~feof(fid))
    line = fgetl(fid);
    linei = linei + 1;
    dblinei = dblinei + 1;
    line_num = str2num(line);
    if line_num(1,2) == ncoeff
        %dblock_i = line_num(1,1);
        datarecord_i = line_num(1,1);
        dblinei = 1;
    end
    if dblinei == 2
        linejd = str2num(line);
        jd1 = linejd(1,1);
        jd2 = linejd(1,2);
        if tjd >= jd1 & tjd <= jd2
            datarecord = [datarecord_i jd1 jd2]';
            % Mercury
            cbv_i = 3;
            decbv(1,1) = linejd(1,3);
        else
            dblinei = -1000000;
        end
    end
    if dblinei > 2 & dblinei < 342
        dbline = str2num(line);
        for i = 1 : length(dbline)
            cbv_coeff = dbline(1,i);
            cbv_i = cbv_i + 1;
            % Mercury
            if (cbv_i >= deformat(1,1)) & (cbv_i < deformat(1,2))
                j = 1;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end
            % Venus
            if (cbv_i >= deformat(1,2)) & (cbv_i < deformat(1,3))
                j = 2;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end
            % Earth-Moon barycenter
            if (cbv_i >= deformat(1,3)) & (cbv_i < deformat(1,4))
                j = 3;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end
            % Mars
            if (cbv_i >= deformat(1,4)) & (cbv_i < deformat(1,5))
                j = 4;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end
            % Jupiter
            if (cbv_i >= deformat(1,5)) & (cbv_i < deformat(1,6))
                j = 5;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end
            % Saturn
            if (cbv_i >= deformat(1,6)) & (cbv_i < deformat(1,7))
                j = 6;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end        
            % Uranus
            if (cbv_i >= deformat(1,7)) & (cbv_i < deformat(1,8))
                j = 7;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end        
            % Neptune
            if (cbv_i >= deformat(1,8)) & (cbv_i < deformat(1,9))
                j = 8;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end
            % Pluto
            if (cbv_i >= deformat(1,9)) & (cbv_i < deformat(1,10))
                j = 9;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end
            % Moon
            if (cbv_i >= deformat(1,10)) & (cbv_i < deformat(1,11))
                j = 10;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end
            % Sun
            if (cbv_i >= deformat(1,11)) & (cbv_i < deformat(1,12))
                j = 11;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end
            % Nutations
            if (cbv_i >= deformat(1,12)) & (cbv_i < deformat(1,13))
                j = 12;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end        
            % Librations
            if (cbv_i == deformat(1,13))
                j = 13;
                decbv_i = cbv_i - (deformat(1,j) - 1);
                decbv(decbv_i,j) = cbv_coeff;
            end
            clear decbv_i
        end
    elseif dblinei == 342
        break
    end
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
