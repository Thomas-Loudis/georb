function [Pnm_norm] = Legendre_functions(phi,nmax) 

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
% Legendre Functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Pnm degree n, order m 
% 
% Computation of the Associated Legendre functions and the Normalized 
% Associated Legendre functions 
% 
% Evaluation of Legendre functions is based on Recurrence relations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Thomas D. Papanikolaou,  AUTH                                   July 2007 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Preallocation 
Pnm_norm = zeros(nmax+1,nmax+1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% computation of the Normalized Associated Legendre functions 
 
%u = sin(phi); 
theta = pi/2 - phi; 
c = cos(theta); 
s = sin(theta); 
 
Poo_norm = 1; 
P10_norm = sqrt(3) * c; 
P11_norm = sqrt(3) * s; 
 
%------------------------------------------------------------------------- 
% Normalized Associated Legendre functions (Pnm_Norm) for order m=0 and 
% variable degree n 
for n = 0: nmax 
    if n==0 
        Pnm_norm(n+1,0+1) = Poo_norm ; 
    elseif n==1 
        Pnm_norm(n+1,0+1) = P10_norm ; 
    else  
        P_f1 = sqrt(2*n-1) * c * Pnm_norm(n-1+1,0+1); 
        Pnm_norm(n+1,0+1) = (sqrt(2*n+1) / n) * ( P_f1 - ( (n-1) / sqrt(2*n-3) ) * Pnm_norm(n-2+1,0+1) ); 
    end 
end 
 
 
% Normalized functions for equal degree n and order m     
for n = 1 : nmax 
    if n == 1 
        Pnm_norm(n+1,n+1) = P11_norm;  
    else    
        Pnm_norm(n+1,n+1) = ( sqrt(2*n+1) / sqrt(2*n) ) * s * Pnm_norm(n-1+1,n-1+1); 
    end 
end 
 
% Normalized functions for variable degree n and order m 
for n = 2 : nmax 
    for m = 1 : n-1 
        % if n >= (m+1) 
            P_f2 = sqrt(2*n-1) * c * Pnm_norm(n-1+1,m+1); 
            P_f3 = sqrt( ( (n-1+m)*(n-1-m) ) / (2*n-3) ) * Pnm_norm(n-2+1,m+1) ; 
            Pnm_norm(n+1,m+1) = ( sqrt(2*n+1) / sqrt( (n+m)*(n-m) ) ) * (P_f2 - P_f3) ; 
        % end 
    end 
end 
