function [RMS, sigma] = rms(dx) 

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
% RMS - root mean square 
% 
% Thomas D. Papanikolaou                                     December 2005 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
[size_dx, size2] = size(dx); 
sum_dx2 = 0; 
for i = 1 : size_dx 
    sum_dx2 = sum_dx2 + dx(i,1)^2 ; 
end 
 
RMS = sqrt(sum_dx2 / size_dx) ; 
 
Xmean = sum_dx2 / size_dx; 
 
sum_dx2 = 0; 
for i = 1 : size_dx 
    sum_dx2 = sum_dx2 + (dx(i,1) - Xmean)^2 ; 
end 
 
sigma = sqrt(sum_dx2 / size_dx); 
 
 
