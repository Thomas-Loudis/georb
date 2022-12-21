function [m1, m2] = wobble_var(mjd,eop,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wobble_var : Wobble variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Wobble variables coputation as function of polar motion variables based
%  on IERS Conventions 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments
% - mjd     : MJD in Terrestrial Time (TT) scale including the fraction of
%             the day 
% - eop:    Earth Orientation Parameters (EOP) data that are required for
%           the orbit arc length
% - dpint:  Number of data points (days) that are required for the EOP
%           interpolation to the computation epoch
%
% Output arguments:
% - m1    : Wobble variable m1 in arcsec
% - m2    : Wobble variable m2 in arcsec
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             8 October 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean Pole coordinates at epoch MJD in mas (milli-arcsec)
[xp_mean_mas, yp_mean_mas] = mean_pole_model(mjd);
xp_mean_arcsec = xp_mean_mas * 10^-3;
yp_mean_arcsec = yp_mean_mas * 10^-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polar Motion coordinates at MJD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eop_mjd:  Array of values of the EOP days refer in seconds (UTC)
eop_mjd = eop(:,1);
xpole   = eop(:,2);
ypole   = eop(:,3);

% Lagrangian Interpolation of EOP parameters is performed at UTC time scale
% Civil date (D,Mh,Yr)
[TT,Dy,Mh,Yr] = MJD_inv(mjd);
% computation of UTC time
[UTC,GPS_time] = time_scales(TT,mjd);
% MJD in UTC time scale
[jd,mjd_int] = MJD_date(UTC,Dy,Mh,Yr);

% Polar motion
[xp_arcsec] = interp_Lagrange(eop_mjd,xpole,mjd_int,dpint);
[yp_arcsec] = interp_Lagrange(eop_mjd,ypole,mjd_int,dpint);

% Conversion from arcsec to radians
%xp = (xp_arcsec / 3600) * (pi / 180);
%yp = (yp_arcsec / 3600) * (pi / 180);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wobble variables (m1, m2) to the polar motion variables (xp, yp)
m1 =    xp_arcsec - xp_mean_arcsec;
m2 = - (yp_arcsec - yp_mean_arcsec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
