function [r,v,chebycoef] = dexxxeph_state(body_ith,tjd,HDfilename,decbv,derecord)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dexxxeph_state :  Computation of state vector based on JPL's DExxx
%                   Planetary and Lunar ephemerides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the state vector for Sun, Moon, Planets based on the
%  JPL's DExxx series ephemerides and the evaluation of Chebychev
%  polynomials approximation.
%
% Input arguments
% - body_ith   : Number of the solar system body according to the DExxx
%                list of the 13 items
% - tjd        : Julian Date number of the required epoch
% - HDfilename : Header file name of the DExxx ephemeris
% - decbv      : DExxx Chebyshev coefficients for the required data record
% - derecord   : DExxx data record's number, start, end JD time
%   >>           derecord = [ datarecord_no JDstart JDend]'
%
% Output arguments:
% - r         : Position vector in Solar Barycentric System   (Km)
% - v         : Velocity vector in Solar Barycentric System   (Km/day)
% - chebycoef : Chebyshev coefficients that were required only
%
% >> r = [ X Y Z]';  v = [Vx Vy Vz]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read DExxx Header file: header.xxx
fid = fopen(HDfilename);
i = 0;
j = 0;
line_no = 0;
linei_group1030 = -1000;
while (~feof(fid))
    line = fgetl(fid);
    i = i + 1;
    if i == 1
        line1num = sscanf(line,'%s %d');
        [sz1 sz2] = size(line1num);
        ncoeff = line1num(sz1,sz2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    linei_group1030 = linei_group1030 + 1;
    if length(line) > 11
        if line(1:12) == 'GROUP   1030'
            linei_group1030 = 1;
        end
    end
    if linei_group1030 == 3
        group1030 = str2num(line);
        tspan = group1030(1,3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
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
% Select required Chebyshev coefficients from data record (array: "decbv")
% and store them into the array "chebycoef". 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DExxx's selected data record for the required epoch tJD
t1 = derecord(2,1);
t2 = derecord(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storing the required Chebychev coefficients from the full data record
if deformat(3,body_ith) == 1
    nchebycoef = deformat(4,body_ith);
    chebycoef = decbv(1:nchebycoef,body_ith);
else
    % Interval of days at computation epoch sicne t1 of the Data record
    tjd_interval = tjd - t1; 
    % Coefficients set of the complete data record period (32 days)
    % Subinterval size of days for the required Celestial body
    subinterval_body = tspan / deformat(3,body_ith);
    % Number of subinterval of the computation epoch
    subinterval_epoch = fix(tjd_interval / subinterval_body) + 1;
    % Read Chebychev coefficients for the computation epoch's subinterval
    chebycoef_indx1 = (subinterval_epoch-1) * deformat(2,body_ith) * 3 + 1;
    chebycoef_indx2 = (subinterval_epoch) * deformat(2,body_ith) * 3;
    chebycoef = decbv(chebycoef_indx1 : chebycoef_indx2,body_ith);
    % JD at the start and end of epoch's subinterval coefficients
    subinterval_t1 = t1 + subinterval_body * (subinterval_epoch-1);
    subinterval_t2 = t1 + subinterval_body * subinterval_epoch;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State vector computation of the selected solar system body in the
% Barycentric reference system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz1 sz2] = size(chebycoef);
chebycoefX = chebycoef(1 : (sz1/3),1);
chebycoefY = chebycoef((sz1/3 + 1) : (2*sz1/3),1);
chebycoefZ = chebycoef((2*sz1/3 + 1) : (3*sz1/3),1);
%[ft,dft] = chebypol(coeff,tjd,subinterval_t1,subinterval_t2)
[X,Vx] = chebypol(chebycoefX,tjd,subinterval_t1,subinterval_t2);
[Y,Vy] = chebypol(chebycoefY,tjd,subinterval_t1,subinterval_t2);
[Z,Vz] = chebypol(chebycoefZ,tjd,subinterval_t1,subinterval_t2);
r = [ X Y Z]';
v = [Vx Vy Vz]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

