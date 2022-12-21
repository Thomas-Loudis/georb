function [eop] = eop_read(filename,mjd)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  eop_read
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read EOP (Earth Orientation Parameters) data series provided by IERS
%  (International Earth Rotation Service and Reference Systems) combined
%  series C04 solution format
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - filename:  EOP C04 data series file name  e.g. 'eopc04_IAU2000.62-now'
% - mjd:       Array of the MJD (Modified Julian Day) numbers of the
%              required dates (MJD in UTC time scale)
%
% Output arguments:
% - eop:      Array of the selected EOP data for the days given by the
%             "mjd" imput argument
%   eop matrix is defined as eop = [mjd x y UT1_UTC dX dY]  nx6 matrix
%   mjd:      MJD refered to 0h in UTC time scale 
%   x,y:      Polar motion coordinates (in seconds) 
%   UT1_UTC:  Difference between UT1 and UTC (in arcsec)
%   dX,dY:    VLBI corrections to the Precession-Nutation model IAU2000
%             (in arcsec)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark: 
%  EOP series that are available by IERS are refered to UTC time scale.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 06/04/2011    Thomas Papanikolaou 
%               Upgrade reading of IERS 08 C04 EOP series format.
% 21/11/2022    Thomas Papanikolaou 
%               Upgrade reading commands and rename to eop_read.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read EOP series data format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nmjd n2] = size(mjd);
clear n2
fid = fopen(filename);
i = 1;
while (~feof(fid))
  line = fgetl(fid);
% FORMAT(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),2(F11.7),2(F12.6))
  param_1 = sscanf(line,'%s %*');
  mjdi = sscanf(line,'%*s%*s%*s %d %*'); 
  if mjdi == mjd(i,1) %1 + i
      mjd_dpi = mjdi;
      x = sscanf(line,'%*s%*s%*s%*s %e %*');
      y = sscanf(line,'%*s%*s%*s%*s%*e %f %*');      
      UT1_UTC = sscanf(line,'%*s%*s%*s%*s%*f%*f %f %*') ;
      dX = sscanf(line,'%*s%*s%*s%*s%*f%*f%*f%*f %f %*') ;
      dY = sscanf(line,'%*s%*s%*s%*s%*f%*f%*f%*f%*f %f %*') ;
      eop(i,:) = [mjd_dpi x y UT1_UTC dX dY];
      i = i + 1;
  end  
  if i > nmjd
      break
  end
end
fclose(fid);
clear i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%