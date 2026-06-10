function [b] = binom_coeff(n,k) 

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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Binomial Coefficients 
% 
% Purpose: 
%   Computation of the binomial coefficients that are reqired in 
%   Adams-Bashforth integration method 
% 
% Input arguments 
% - n,k:  Elements 
% 
% Output arguments 
% - b:   Binomial coefficient of elements n,k according to the mathematical 
%        scheme b_nk = (n k), 
%        where n and k stand for the up and down elements respectively. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Thomas D. Papanikolaou                                         May 2010 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
factor = 1; 
for i = 0 : k-1 
    factor = factor * (n - i); 
end 
 
b = factor / factorial(k); 
