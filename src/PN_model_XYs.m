function [XYs] = PN_model_XYs(mjd, IAU_model)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  PN_model_XYs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  IAU Precession-Nutation model : X,Y,s parameters
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - IAU_model: IAU Precession-Nutation model: IAU2000A or IAU2006/2000A: 
%              Variable values : 2000 or 2006
% - mjd:       Array of the MJD (Modified Julian Day) numbers of the
%              required dates (MJD in UTC time scale)
% - filename:  EOP C04 data series file name  e.g. 'eopc04_IAU2000.62-now'
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
% Thomas Loudis Papanikolaou                               22 November 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% XYs matrix format :   MJD  X  Y  s
test = strcmp(IAU_model,'IAU2000A');
if test == 1
    [XYs_model_array] = pn_iau2000A(mjd);
end
test = strcmp(IAU_model,'IAU2006');
if test == 1
    [XYs_model_array] = pn_iau2006(mjd);  
end

[n m] = size(XYs_model_array);
[sz1 sz2] = size(mjd);
XYs = zeros(sz1,4);
xys_indx = 1;
for i_XYs = 1 : n    
    XYZ_epoch = XYs_model_array(i_XYs,:);
    mjd_PN = XYZ_epoch(1,1);
    mjd_in = fix(mjd(1,1));
    if abs(mjd_in - mjd_PN) < 10^-8     
        for i = 1 : sz1
            mjd_epoch = fix(mjd(i,1));
            XYZ_epoch = XYs_model_array(i_XYs,:);
            XYs(xys_indx,:) = [mjd_epoch XYZ_epoch(1,2) XYZ_epoch(1,3) XYZ_epoch(1,4)];            
            xys_indx = xys_indx + 1;
            i_XYs = i_XYs + 1;
        end
        break
    end    
end

