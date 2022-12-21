function [Cnm_tides,Snm_tides] = tides_add2(Cnm,Snm,dCnm,dSnm,n_max)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tides_add2 :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Add tide corrections to the spherical harmonics coefficients
%  i.e. add dCnm,dSnm corrections to the Cnm Snm arrays
%
% Input arguments
% - Cnm   :  Gravity model C spherical harmonic coefficients
% - Snm   :  Gravity model S spherical harmonic coefficients
% - dCnm  :  Cnm coefficietns tide corrections
% - dSnm  :  Snm coefficietns tide corrections
% - n_max :  Degree/Order limit for the tide corrections to be added 
%
% Output arguments:
% - Cnm_tides : Cnm matrix plus Tides corrections
% - Snm_tides : Snm matrix plus Tides corrections
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
%  Computed dCnm and dSnm are stored into lower triangular matrices.
%  Coefficient dCnm corresponds to matrix element dCnm(n+1,m+1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                        June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 06/06/2012  Function's upgrade and renamed to tides_add2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Cnm maximum degree
[n1 n2] = size(Cnm);
Nmax_Cnm = n1 - 1;
clear n1 n2

% dCnm maximum degree
[n1 n2] = size(dCnm);
Nmax_dCnm = n1 - 1;
clear n1 n2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum Degree comparison between Cnm & dCnm
if Nmax_Cnm < Nmax_dCnm
    Nmax = Nmax_Cnm;
else
    Nmax = Nmax_dCnm;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Degree/Order limit from input argument
if n_max > 1
   Nmax = n_max;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation
Cnm_tides = zeros(Nmax+1,Nmax+1);
Snm_tides = zeros(Nmax+1,Nmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add coefficients corrections
for n_i = 2 : Nmax
    for m_i = 0 : n_i
        Cnm(n_i + 1,m_i + 1) = Cnm(n_i + 1,m_i + 1) + dCnm(n_i + 1,m_i + 1);
        Snm(n_i + 1,m_i + 1) = Snm(n_i + 1,m_i + 1) + dSnm(n_i + 1,m_i + 1);
    end
end
clear Nmax n_i m_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Final matrices
Cnm_tides = Cnm;
Snm_tides = Snm;
