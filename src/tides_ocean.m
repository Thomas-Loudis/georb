function [dCnm,dSnm] = tides_ocean(n_max,m_max,mjd,eop,dpint,DelaunayNf,dCnm_plus,dSnm_plus,dCnm_minus,dSnm_minus, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tides_ocean : Ocean Tides perturbation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Ocean Tides perturbations computation according to the formulas of the
%  IERS Conventions 2010. This function implements the computation of the
%  variations of the Stokes coefficients based on the ocean tide model
%  FES2004.
%
% Input arguments
% - n_max      : Degree limit of the ocean tides coefficients
% - m_max      : Order limit of the ocean tides coefficients
% - mjd          : MJD in Terrestrial Time (TT) scale including the fracti-
%                  on of the day 
% - eop          : Earth Orientation Parameters (EOP) data that are
%                  required for the orbit arc length
% - dpint        : Number of data points (days) that are required for the 
%                  EOP interpolation to the computation epoch
% - DelaunayNf : Delaunay variables multipliers for FES2004 tidal waves
% - dCnm_plus  : FES2004 coefficients to be used in Eq.6... IERS Conv.2010
%   dSnm_plus    Coefficients are stored in these four 3D arrays
%   dCnm_minus
%   dSnm_minus
%
% Output arguments:
% - dCnm  :  Gravity field's Cnm corrections matrix due to Ocean Tides
% - dSnm  :  Gravity field's Snm corrections matrix due to Ocean Tides
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
%  Computed dCnm and dSnm are stored into lower triangular matrices.
%  Coefficient dCnm corresponds to matrix element dCnm(n+1,m+1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
%  23/06/2012   Upgraded for compatibility with modification of function
%               tides_fes2004.m
%               Computation time has now been reduced dramatically for
%               ocean tides perturbations.
% 
% 01/02/2019  Dr. Thomas Papanikolaou
%             Correction to the basic Equation of the ocean tides (IERS Conv. 2010 Eq. 6.15)
%             Correction to the multiplier of the theta_f angle computation 
% 
% 15/01/2025  Dr. Thomas Loudis Papanikolaou
%             Code optimization (CPU time minimization) 
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocation
% dCnm = zeros(n_max+1,n_max+1);
% dSnm = zeros(n_max+1,n_max+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delaunay variables (in radians)
[F1,F2,F3,F4,F5] = delaunay_variables(mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Greenwich Mean Sidereal Time (GMST) in radians
[thetag] = iers_gmst(mjd,eop,dpint, orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thetaf computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delaunay and Doodson multipliers :: DelaunayNf matrix 
% DelaunayNf = [Doodson_Number N1 N2 N3 N4 N5 n1 n2 n3 n4 n5 n6]
[Nfrq sz2] = size(DelaunayNf);
thetaf_matrix = zeros(Nfrq,1);
thetaf_vec = zeros(1,Nfrq);
thetaf_3d = zeros(1,1,Nfrq);
cos_thetaf_3d = zeros(1,1,Nfrq);
sin_thetaf_3d = zeros(1,1,Nfrq);
for ifrq = 1 : Nfrq
    %thetaf = Doodson_n1 * (thetag + pi) - DelaunayNf(ifrq,2:6) * [F1 F2 F3 F4 F5]';
    % thetaf (in radians)
    
    % Delaunay_freq = DelaunayNf(ifrq,1);
    Delaunay_N1   = DelaunayNf(ifrq,2);
    Delaunay_N2   = DelaunayNf(ifrq,3);
    Delaunay_N3   = DelaunayNf(ifrq,4);
    Delaunay_N4   = DelaunayNf(ifrq,5);
    Delaunay_N5   = DelaunayNf(ifrq,6);
    Doodson_number = DelaunayNf(ifrq,1);
    Doodson_n1     = DelaunayNf(ifrq,7);
    
    %DelaunayNf(ifrq,2:6) * [F1 F2 F3 F4 F5]'
    Delaunay_sum = Delaunay_N1 * F1 + ...
        Delaunay_N2 * F2 + ...
        Delaunay_N3 * F3 + ...
        Delaunay_N4 * F4 + ...
        Delaunay_N5 * F5 ;
    thetaf = Doodson_n1 * (thetag + pi) - Delaunay_sum;
    thetaf_matrix(ifrq,1) = thetaf;
    thetaf_vec(1,ifrq) = thetaf;
    thetaf_3d(1,1,ifrq) = thetaf;
    cos_thetaf_3d(1,1,ifrq) = cos(thetaf);
    sin_thetaf_3d(1,1,ifrq) = sin(thetaf);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dCnm_fi = (dCnm_plus + dCnm_minus) .* cos_thetaf_3d + (dSnm_plus + dSnm_minus) .* sin_thetaf_3d;
dSnm_fi = (dSnm_plus - dSnm_minus) .* cos_thetaf_3d - (dCnm_plus - dCnm_minus) .* sin_thetaf_3d;

dCnm_f = sum(dCnm_fi,3);
dSnm_f = sum(dSnm_fi,3);
% IERS Conventions update 10/08/2012 Section 6.3.2
% m = 0 :: dSnm = 0;
[d1, d2] = size(dSnm_f);
dSnm_f(:,0+1) = zeros(d1,1);

% Preallocation
% dCnm = zeros(n_max+1,n_max+1);
% dSnm = zeros(n_max+1,n_max+1);

% Final Matrices
dCnm = dCnm_f(1 : n_max+1 , 1 : n_max+1);
dSnm = dSnm_f(1 : n_max+1 , 1 : n_max+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
