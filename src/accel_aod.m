function [fx,fy,fz] = accel_aod(r,n_max,m_max,GM,ae, Cnm, Snm, MJD)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleration due to Atmosphere and Ocean De-Aliasing (AOD) effects data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Computation of the acceleration in space based on the processing of the
%   AOD data (Stokes coefficients) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - r           :   position vector in Terrestrial Reference System (ITRS)
%                   r = [x y z]'
% - n_max       :   Maximum degree of the spherical harmonics expansions series
% - m_max       :   Maximum order of the spherical harmonics expansions series
% - Cnm         :   AOD C Stokes coefficients, 3-dimensional matrix for data sets of the day (data set per 3 hours)   
% - Snm         :   AOD S Stokes coefficients, 3-dimensional matrix for data sets of the day (data set per 3 hours)
% - mjd         :   Modified Julian Day number including fraction of the day
%
% Output arguments:
% - fx,fy,fz    :   Acceleration's cartesian components in ITRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             2 October 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
%  9/12/2022    Thomas Loudis Papanikolaou
%               Code upgrade and call of the new function harmonics_partials1.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seconds since 00h in TT
TT = (MJD - fix(MJD)) * 86400;
% MJD in UTC time 
[UTC, GPS_time] = time_scales(TT,MJD);
MJD_utc = MJD + (UTC - TT) / 86400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stokes coefficients interpolation at current epoch (linear interpolation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iput Epoch in hours since 00h
MJD_hours = (MJD - fix(MJD)) * 24;

% Periods of 3 hours at input epoch
MJD_3h = MJD_hours / 3;

% AOD: data sets' Numbers for interpolation
AOD_t0_dataset = fix(MJD_3h) + 1;
AOD_t1_dataset = AOD_t0_dataset + 1;
if AOD_t1_dataset > 8
AOD_t0_dataset = AOD_t0_dataset - 1;
AOD_t1_dataset = AOD_t1_dataset - 1;    
end

% AOD data sets' Start and End time in hours since 00h
AOD_t0_h = fix(MJD_3h) * 3;
AOD_t1_h = AOD_t0_h + 3;

% Linear Interpolation of Stokes coefficients

% dt between two data sets in hours
delta_t_datasets = 3;

% dt of input epoch in hours since interplation start epoch
delta_t_epoch = MJD_hours - AOD_t0_h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients interpolation
Cnm_mjd = zeros(n_max+1,n_max+1);
Snm_mjd = zeros(n_max+1,n_max+1);
for n = 0 : n_max
for m = 0 : n
% C coefficients
Cto = Cnm(n+1, m+1, AOD_t0_dataset);
Ct1 = Cnm(n+1, m+1, AOD_t1_dataset);
% Cnm interpolated value
Cnm_mjd(n+1,m+1) = Cto + ( (Ct1 - Cto) / delta_t_datasets ) * delta_t_epoch ;

% S coefficients
Sto = Snm(n+1, m+1, AOD_t0_dataset);
St1 = Snm(n+1, m+1, AOD_t1_dataset);
% Snm interpolated value
Snm_mjd(n+1,m+1) = Sto + ( (St1 - Sto) / delta_t_datasets ) * delta_t_epoch ;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Harmonic series partilas w.r.t. position vector (spherical and
% cartesian coordiantes)
[partials_rpl, partials_xyz] = potential_partials_1st(r,n_max,m_max,GM,ae,Cnm_mjd,Snm_mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cartesian counterparts of the acceleration vector 
fx = partials_xyz(1,1);
fy = partials_xyz(2,1);
fz = partials_xyz(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
