function [dCnm,dSnm] = tides_ocean(n_max,m_max,mjd,eop,dpint,DelaunayNf,dCnm_plus,dSnm_plus,dCnm_minus,dSnm_minus)


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% dCnm = zeros(n_max,n_max,Nfrq);
% dSnm = zeros(n_max,n_max,Nfrq);
% Preallocation
dCnm = zeros(n_max+1,n_max+1);
dSnm = zeros(n_max+1,n_max+1);

% loop_before = 1
% [N_dCnm_plus, m_dCnm_plusm, N_freq] = size(dCnm_plus)
% [N_dCnm, m_dCnm] = size(dCnm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delaunay variables (in radians)
[F1,F2,F3,F4,F5] = delaunay_variables(mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Greenwich Mean Sidereal Time (GMST) in radians
[thetag] = iers_gmst(mjd,eop,dpint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thetaf computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delaunay and Doodson multipliers :: DelaunayNf matrix 
% DelaunayNf = [Doodson_Number N1 N2 N3 N4 N5 n1 n2 n3 n4 n5 n6]
[Nfrq sz2] = size(DelaunayNf);
thetaf_matrix = zeros(Nfrq,1);
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Nfrq sz2] = size(DelaunayNf);
clear sz2
for n = 2 : n_max
    if n > m_max
        m_limit = m_max;
    else
        m_limit = n;
    end
    for m = 0 : m_limit
        dCnm_f = 0;
        dSnm_f = 0;
        for ifrq = 1 : Nfrq  
%             % thetaf (in radians)
%             %thetaf = Doodson_n1 * (thetag + pi) - DelaunayNf(ifrq,2:6) * [F1 F2 F3 F4 F5]';
%             % Delaunay_freq = DelaunayNf(ifrq,1);
%             Delaunay_N1   = DelaunayNf(ifrq,2);
%             Delaunay_N2   = DelaunayNf(ifrq,3);
%             Delaunay_N3   = DelaunayNf(ifrq,4);
%             Delaunay_N4   = DelaunayNf(ifrq,5);
%             Delaunay_N5   = DelaunayNf(ifrq,6);
%             Doodson_number = DelaunayNf(ifrq,1); 
%             Doodson_n1     = DelaunayNf(ifrq,7);            
%             %DelaunayNf(ifrq,2:6) * [F1 F2 F3 F4 F5]'
%             Delaunay_sum = Delaunay_N1 * F1 + ... 
%                            Delaunay_N2 * F2 + ...
%                            Delaunay_N3 * F3 + ...
%                            Delaunay_N4 * F4 + ...
%                            Delaunay_N5 * F5 ;                           
%             thetaf = Doodson_n1 * (thetag + pi) - Delaunay_sum;             
            thetaf = thetaf_matrix(ifrq,1);
            dCnm_fi = (dCnm_plus(n+1,m+1,ifrq) + dCnm_minus(n+1,m+1,ifrq)) * cos(thetaf) + (dSnm_plus(n+1,m+1,ifrq) + dSnm_minus(n+1,m+1,ifrq)) * sin(thetaf);
            dSnm_fi = (dSnm_plus(n+1,m+1,ifrq) - dSnm_minus(n+1,m+1,ifrq)) * cos(thetaf) - (dCnm_plus(n+1,m+1,ifrq) - dCnm_minus(n+1,m+1,ifrq)) * sin(thetaf);
            dCnm_f = dCnm_f + dCnm_fi;
            dSnm_f = dSnm_f + dSnm_fi;
        end
        % IERS Conventions update 10/08/2012 Section 6.3.2
        if m == 0
          dSnm_f = 0;
        end
        dCnm(n+1,m+1) = dCnm_f;
        dSnm(n+1,m+1) = dSnm_f;
        dCnm_f = 0;
        dSnm_f = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop_after = 1
% [N_dCnm_plus, m_dCnm_plusm, N_freq] = size(dCnm_plus)
% [N_dCnm, m_dCnm] = size(dCnm)
