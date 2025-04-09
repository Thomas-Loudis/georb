function [dCnm,dSnm] = tides_pole_ocean(mjd,eop,dpint,GM_Earth,radius_Earth,n_max,m_max, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tides_pole_ocean : Ocean Pole Tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Ocean Pole Tides based on self-consistent equilibrium model as presented
%  by Desai (2002) and adopted by IERS Conventions 2010 update.  
%
% Input arguments
% - mjd     : MJD in Terrestrial Time (TT) scale including the fraction of
%             the day 
% - eop:    Earth Orientation Parameters (EOP) data that are required for
%           the orbit arc length
% - dpint:  Number of data points (days) that are required for the EOP
%           interpolation to the computation epoch
%
% Output arguments:
% - dCnm    : Cnm corrections matrix
% - dSnm    : Snm corrections matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
%  Computed dCnm and dSnm are stored into lower triangular matrices.
%  Coefficient dCnm corresponds to matrix element dCnm(n+1,m+1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             8 October 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 29/08/2024, Thomas Loudis Papanikolaou
%             Upgrade to include higher degree terms based on Desai (2002)
%             model and coeeficients data
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Forces model structure matrix 
oceanpoletide_struct = orbit_model_struct.oceanpoletide;
loadlovenumbers_struct = orbit_model_struct.loadlovenumbers;
TAI_UTC_table = orbit_model_struct.TAI_UTC_table;

% dCnm = zeros(3,3);
% dSnm = zeros(3,3);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Wobble variables (m1, m2) in arcsec 
% [m1, m2] = wobble_var(mjd,eop,dpint);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Solid Earth Pole Tide: changes to geopotential coefficients C21 and S21
% dC21 = -2.1778 * 10^-10 * (m1 - 0.01724 * m2);
% dS21 = -1.7232 * 10^-10 * (m2 - 0.03365 * m1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % dCnm and dSnm corresponding matrices
% dCnm(2+1,1+1) = dC21;
% dSnm(2+1,1+1) = dS21;


% Preallocation
dCnm = zeros(n_max+1, n_max+1);
dSnm = zeros(n_max+1, n_max+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wobble variables (m1, m2) in arcsec 
[m1, m2] = wobble_var(mjd,eop,dpint,TAI_UTC_table);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Love numbers (LLN) kn' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kn_lln = orbit_model_forces_glob.loadlovenumbers_kn;
kn_lln = loadlovenumbers_struct.kn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GM = GM_Earth;
ae = radius_Earth;

% G gravitational constant
G_earth = 6.6730000000e-11; % m^3/(kg*s^2)

% Earth's  mean angular velocity
Omega_earth = 7.2921150000e-05;  % rad/s

% Earth's mean gravity on surface
ge = 9.7803278000e+00;  % m/s^2

% water density Kg/m^3
rho_w = 1025;

% gama2 real part
gama2_R = 0.6870;
% gama2 imaginary part
gama2_I = 0.0036;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficient for Conversion from arcsec to radians
arcsec2rad = pi / (180 * 3600);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rn coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant part of Rn coefficients
Rn_term0 = ( (Omega_earth^2 * ae^4) / GM ) * ( (4*pi * G_earth * rho_w) / ge );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients Anm and Bnm including real and imaginary part based on Desai
% as provided by IERS data file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AnmR real part
% AnmI Imaginary part
Anm_R = oceanpoletide_struct.desaioplecoef_Anm_R;
Anm_I = oceanpoletide_struct.desaioplecoef_Anm_I ;
Bnm_R = oceanpoletide_struct.desaioplecoef_Bnm_R;
Bnm_I = oceanpoletide_struct.desaioplecoef_Bnm_I;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IERS Conventions 2010 (update 10/08/2012), Section 6.5

% Truncated matrices
Anm_R_trunc = Anm_R(1 : n_max+1 , 1 : n_max+1);
Anm_I_trunc = Anm_I(1 : n_max+1 , 1 : n_max+1);
Bnm_R_trunc = Bnm_R(1 : n_max+1 , 1 : n_max+1);
Bnm_I_trunc = Bnm_I(1 : n_max+1 , 1 : n_max+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rn coefficients
Rn_matrix = zeros(n_max+1,1);
for n = 2 : n_max
    Rn_matrix(n+1,1) = Rn_term0 * ( (1 + kn_lln(n+1) ) / (2*n+1) );        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coefficients
coef_R = (m1 * arcsec2rad * gama2_R + m2 * arcsec2rad * gama2_I);
coef_I = (m2 * arcsec2rad * gama2_R - m1 * arcsec2rad * gama2_I);

% sum of matrices
term_A = ( Anm_R_trunc * coef_R + Anm_I_trunc * coef_I );
term_B = ( Bnm_R_trunc * coef_R + Bnm_I_trunc * coef_I );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Harmonics matrices
dCnm = term_A .* Rn_matrix;
dSnm = term_B .* Rn_matrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
