function [Cnm_sum,Snm_sum] = harmonics_sum(Cnm,Snm,dCnm,dSnm,n_max)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shc_add :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Add harmonics coefficients corrections 
%
% Input arguments
% - Cnm   :  Gravity model C spherical harmonic coefficients
% - Snm   :  Gravity model S spherical harmonic coefficients
% - dCnm  :  Cnm coefficiets corrections
% - dSnm  :  Snm coefficiets corrections
% - n_max :  Degree/Order limit for the corrections to be added 
%
% Output arguments:
% - Cnm_sum : Cnm matrix plus corrections
% - Snm_sum : Snm matrix plus corrections
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                        June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 06/06/2012  Function's upgrade and renamed to tides_add2.m
% 30/03/2023  Function revision and rename to be used in gravity parameter
%             estimation
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
Cnm_sum = zeros(Nmax+1,Nmax+1);
Snm_sum = zeros(Nmax+1,Nmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add coefficients corrections
for n_i = 0 : Nmax
    for m_i = 0 : n_i
        Cnm(n_i + 1,m_i + 1) = Cnm(n_i + 1,m_i + 1) + dCnm(n_i + 1,m_i + 1);
        Snm(n_i + 1,m_i + 1) = Snm(n_i + 1,m_i + 1) + dSnm(n_i + 1,m_i + 1);
    end
end
clear Nmax n_i m_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Final matrices
Cnm_sum = Cnm;
Snm_sum = Snm;
