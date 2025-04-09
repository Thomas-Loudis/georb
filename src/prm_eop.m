function [EOP_data_array, IAU_PN_XYs_matrix] = prm_eop(eopfilename, eop_interp_points, time_period, MJDo, IAU_PN_model, TAI_UTC_table) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prm_eop : Read Earth Orientation Parameters data and select the
% parameters for the period of the time series analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - eopfilename        : IERS EOP data file name
% - eop_interp_points  : Number of data points (days) required for the EOP
%                        interpolation to the computation epoch
% - time_period : Overall period of time analysis in Seconds
%
% Output arguments:
% - EOP_data_array     : Earth Orientation Parameters (EOP) data array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Papanikolaou                                         January  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 28/05/2022   Dr. Thomas Loudis Papanikolaou  
%              prm3 function has been modified and renamed to prm_eop 
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numebr of data points (days) used for interpolation of EOP data
dpint = eop_interp_points;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJD of initial epoch's date in UTC time scale
[TT,day,month,year] = MJD_inv(MJDo);
[UTC,GPS_time] = time_scales(TT,MJDo,TAI_UTC_table);
% dt UTC-TT in days
dt_UTC_TT = (UTC - TT) / (24*3600);
mjdo_UTC = MJDo + dt_UTC_TT;
% MJD number of the initial epoch's day (Remark: converse MJD to integer)
mjdo_UTC = fix(mjdo_UTC);
if mjdo_UTC < fix(MJDo)
    mjdo_UTC = fix(MJDo);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time period length in Days
arcd = fix(time_period / (24 * 60 * 60)) + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Half of dpint
dpint2 = dpint/2;
% Integer of dpint2
dpint2_intg = fix(dpint2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Array of MJD numbers of EOP data points that are required for the orbit
% arc length: mjd_dp  dpx1 matrix 
if dpint2_intg - dpint2 < 10^-8
    mjd1 = mjdo_UTC - (dpint2 - 1);
    dp = arcd + (dpint2-1) + (dpint2);
    mjd_dp = zeros(dp,1);
    mjd_dp(1,1) = mjd1;
    for i = 2 : dp
        mjd_dp(i,1) = mjd1 + i - 1;
    end
elseif dpint2_intg - dpint2 ~= 0    
    mjd1 = mjdo_UTC - dpint2_intg;
    dp = arcd + dpint2_intg + dpint2_intg;
    mjd_dp = zeros(dp,1);
    mjd_dp(1,1) = mjd1;
    for i = 2 : dp
        mjd_dp(i,1) = mjd1 + i - 1;
    end        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read EOP series data file
[eopdat] = eop_read(eopfilename,mjd_dp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read IAU2000A precession-nutation model :: X,Y,s
eopdat_mjd = eopdat(:,1);
[XYs] = PN_model_XYs(eopdat_mjd, IAU_PN_model);
XYs_IAU200A = XYs;
IAU_PN_XYs_matrix = XYs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EOP_data_array = eopdat; 
