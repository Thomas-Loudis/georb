function [fx,fy,fz] = accel_aod(r,n_max,m_max,GM,ae, Cnm, Snm, MJD, legendre_functions_struct, aod_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleration of Atmosphere and Ocean De-Aliasing (AOD) effects
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
% 16/03/2025    Thomas Loudis Papanikolaou
%               Code upgrade  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Harmonics Synthesis start degree
n_min = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First and Last epochs
mjd_first_epoch_TT = aod_struct.mjd_first_epoch_TT;
mjd_last_epoch_TT  = aod_struct.mjd_last_epoch_TT;
data_sets_number = aod_struct.data_sets_number;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Periods of 3 hours at input epoch since datasets first epoch
MJD_3h = (MJD - mjd_first_epoch_TT) * 24 / 3;

% AOD: Find data sets' Numbers to be used for interpolation
[d1,d2,d3] = size(Cnm);
AOD_t0_dataset = fix(MJD_3h) + 1;
AOD_t1_dataset = AOD_t0_dataset + 1;
if AOD_t1_dataset > d3    
AOD_t0_dataset = d3 - 1;
AOD_t1_dataset = d3    ;
end

% Data Sets start and end epochs in hours since data first epoch
AOD_t0_h = (AOD_t0_dataset - 1) * 3;
AOD_t1_h = AOD_t0_h + 3;

% Data Sets start and end epochs in MJD
MJD_AOD_to = AOD_t0_h / 24 + (mjd_first_epoch_TT);
MJD_AOD_t1 = AOD_t1_h / 24 + (mjd_first_epoch_TT);

% dt between two data sets in hours
delta_t_datasets = 3;

% dt of input epoch in hours since interpolation start epoch
delta_t_epoch = (MJD - MJD_AOD_to) * 24;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stokes coefficients interpolation at current epoch (linear interpolation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interp_apporach = 1; 
if interp_apporach == 1
% 1. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C coefficients
Cto = Cnm(:,:, AOD_t0_dataset);
Ct1 = Cnm(:,:, AOD_t1_dataset);
% Cnm interpolated value
Cnm_mjd = Cto + (Ct1 - Cto) * (delta_t_epoch / delta_t_datasets ) ;

% S coefficients
Sto = Snm(:,:, AOD_t0_dataset);
St1 = Snm(:,:, AOD_t1_dataset);
% Snm interpolated value
Snm_mjd = Sto + (St1 - Sto) * ( delta_t_epoch / delta_t_datasets ) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif interp_apporach == 2
% 2. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients interpolation (2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ti = MJD;
t1 = MJD_AOD_to;
t2 = MJD_AOD_t1;

delta_tit2 = t2 - ti;
delta_t1ti = ti - t1;
delta_t1t2 = t2 - t1;

% C coefficients
Cto = Cnm(:,:, AOD_t0_dataset);
Ct1 = Cnm(:,:, AOD_t1_dataset);
% Cnm interpolated value
Cnm_mjd = Cto * ( delta_tit2 / delta_t1t2) + Ct1 * ( delta_t1ti / delta_t1t2);

% S coefficients
Sto = Snm(:,:, AOD_t0_dataset);
St1 = Snm(:,:, AOD_t1_dataset);
% Snm interpolated value
Snm_mjd = Sto * ( delta_tit2 / delta_t1t2) + St1 * ( delta_t1ti / delta_t1t2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Harmonic series partilas w.r.t. position vector (spherical and
% cartesian coordiantes)
[partials_rpl, partials_xyz] = potential_partials_1st(r,n_max,m_max,GM,ae,Cnm_mjd,Snm_mjd, legendre_functions_struct, n_min);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cartesian counterparts of the acceleration vector 
fx = partials_xyz(1,1);
fy = partials_xyz(2,1);
fz = partials_xyz(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
