function [E_anomaly] = kepler_eq(M,ecc) 

%-------------------------------------------------------------------------- 
% Copyright (C) 2007-present  Thomas (Loudis) Papanikolaou  
%  
% GEORB :: 'Gravity and Precise Orbit Determination system' 
%  
% GEORB is a software for Precise Orbit Determination (POD), gravity field  
% recovery and mission design 
%  
% GEORB is licensed under the GNU Affero General Public License v3.0 while  
% for information regarding dual-license options visit the official website 
% www.georb.gr  
%-------------------------------------------------------------------------- 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Kepler's Equation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%   Solution of Kepler's equation that is based on Newton method 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - M: Mean amomaly in degrees 
% - e: eccentricity 
% 
% Output arguments: 
% - E: Eccentric anomaly in degrees 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Thomas D. Papanikolaou                                     November 2007 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified 
% 8/11/2022    Dr. Thomas Papanikolaou 
%              Code minor upgarde  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
% converse in radians 
M = M * (pi/180) ; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Initial approximate value for Eccentric anomaly (E) 
%  Reccomendation: 
%  - for small eccentricity: Eo = M 
%  - for high eccentricity (e>0.8): Eo = pi 
if ecc < 0.8     
    if M==0 
        Eo = M + 0.001; 
    else 
        Eo = M; 
    end     
else     
    Eo = pi;     
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% number of iterations 
iter = 7; 
% E = zeros(1,iter+1) 
for i = 1 : iter 
    if i == 1 
        E(i) = Eo; 
    end     
    % auxiliary function 
    f_E = E(i) - ecc * sin( E(i) ) - M;  
 
    E(i+1) = E(i) - f_E / ( 1 - ecc * cos(E(i)) ) ; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
[n1 n2] = size(E); 
% Converse in degrees 
E_anomaly = E(n1,n2) * (180 / pi); 
