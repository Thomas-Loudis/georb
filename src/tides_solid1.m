function [dCnm,dSnm] = tides_solid1(GMearth,Re,GMmoon,rmoon,GMsun,rsun)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tides_solid1 : Solid Earth Tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Solid Earth Tides based on IERS Conventions 2010. This function impleme-
%  nts the Step1 according to IERS Conventions 2010 which refer to the com-
%  putation of the frequency indepedent corrections to geopetential's
%  spherical harmonic coefficients. 
%
% Input arguments
% - GMearth : Earth gravity constant  (m^3/sec^2)
% - Re      : radius  (meters)
% - GMmoon  : Moon gravity constant  (m^3/sec^2)
% - rmoon   : Moon body-fixed geocentric position vector (meters)
% - GMsun   : Sun gravity constant  (m^3/sec^2)
% - rsun    : Sun body-fixed geocentric position vector (meters)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dCnm = zeros(5,5);
dSnm = zeros(5,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Knm matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elastic Earth
knm(2+1,0+1) = 0.29525;
knm(2+1,1+1) = 0.29470;
knm(2+1,2+1) = 0.29801;
knm(3+1,0+1) = 0.093;
knm(3+1,1+1) = 0.093;
knm(3+1,2+1) = 0.093;
knm(3+1,3+1) = 0.094;
% Knm(+) matrix
knm_plus(2+1,0+1) = -0.00087;
knm_plus(2+1,1+1) = -0.00079;
knm_plus(2+1,2+1) = -0.00057;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Anelastic Earth
%knm(2+1,0+1) = 0.30190;
%knm(2+1,1+1) = 0.29830;
%knm(2+1,2+1) = 0.30102;
% Knm(+) matrix
%knm_plus(2+1,0+1) = -0.00089;
%knm_plus(2+1,1+1) = -0.00080;
%knm_plus(2+1,2+1) = -0.00057;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of spherical coordinates in radians and distance in meters
% computation of normalized associated Legendre functions
% Moon
[lamda_moon,phi_moon,lmoon] = lamda_phi(rmoon);
[Pnm_norm_moon] = Legendre_functions(phi_moon,3);
% Sun
[lamda_sun,phi_sun,lsun] = lamda_phi(rsun);
[Pnm_norm_sun] = Legendre_functions(phi_sun,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of dCnm and dSnm terms for degree and order 3
for n = 2 : 3
    for m = 0 : n
        % dCnm term
        dCnm_moon = (GMmoon / GMearth) * (Re / lmoon)^(n+1) * Pnm_norm_moon(n+1,m+1) * cos(m*lamda_moon);
        dCnm_sun  = (GMsun / GMearth) * (Re / lsun)^(n+1) * Pnm_norm_sun(n+1,m+1) * cos(m*lamda_sun);
        dCnm(n+1,m+1) = ( knm(n+1,m+1) / (2*n + 1) ) * (dCnm_moon + dCnm_sun);
        % dSnm term
        dSnm_moon = (GMmoon / GMearth) * (Re / lmoon)^(n+1) * Pnm_norm_moon(n+1,m+1) * sin(m*lamda_moon);
        dSnm_sun  = (GMsun / GMearth) * (Re / lsun)^(n+1) * Pnm_norm_sun(n+1,m+1) * sin(m*lamda_sun);
        dSnm(n+1,m+1) = ( knm(n+1,m+1) / (2*n + 1) ) * (dSnm_moon + dSnm_sun);
    end
end
%clear n m
n = 0; m= 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of dCnm and dSnm terms for degree n=4 and order m=0,1,2
for m = 0 : 2
    % dCnm terms
    dCnm_moon = (GMmoon / GMearth) * (Re / lmoon)^3 * Pnm_norm_moon(2+1,m+1) * cos(m*lamda_moon);
    dCnm_sun  = (GMsun / GMearth) * (Re / lsun)^3 * Pnm_norm_sun(2+1,m+1) * cos(m*lamda_sun);
    dCnm(4+1,m+1) = ( knm_plus(2+1,m+1) / 5 ) * (dCnm_moon + dCnm_sun);
    % dSnm terms
    dSnm_moon = (GMmoon / GMearth) * (Re / lmoon)^3 * Pnm_norm_moon(2+1,m+1) * sin(m*lamda_moon);
    dSnm_sun  = (GMsun / GMearth) * (Re / lsun)^3 * Pnm_norm_sun(2+1,m+1) * sin(m*lamda_sun);
    dSnm(4+1,m+1) = ( knm_plus(2+1,m+1) / 5 ) * (dSnm_moon + dSnm_sun);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
