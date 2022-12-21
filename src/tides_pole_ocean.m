function [dCnm,dSnm] = tides_pole_ocean(mjd,eop,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tides_pole_solidearth : Solid Earth Pole Tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Solid Earth Pole Tides based on IERS Conventions 2010. 
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


dCnm = zeros(3,3);
dSnm = zeros(3,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wobble variables (m1, m2) in arcsec 
[m1, m2] = wobble_var(mjd,eop,dpint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solid Earth Pole Tide: changes to geopotential coefficients C21 and S21
dC21 = -2.1778 * 10^-10 * (m1 - 0.01724 * m2);
dS21 = -1.7232 * 10^-10 * (m2 - 0.03365 * m1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dCnm and dSnm corresponding matrices
dCnm(2+1,1+1) = dC21;
dSnm(2+1,1+1) = dS21;
