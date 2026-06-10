function [Rrtn,er,et,en] = orbital_transf(r_crf,v_crf) 

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
% Transformation matrix from inertial to orbital frame 
%  
% Purpose: 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - r_crf: position vector in the inertila frame r_crf = [x y z]' 
%  
% Output arguments: 
% - Rrtn : Transformation matrix between the inertial and orbital frames 
%          radial, along-track and cross-track components 
% Remark: transformation matrix is applied as follows 
% [ar at an]' = Rrtn * [ax ay az]' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas D. Papanikolaou                                    August 2013 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
r_length = sqrt(r_crf(1,1)^2 + r_crf(2,1)^2 + r_crf(3,1)^2); 
v_length = sqrt(v_crf(1,1)^2 + v_crf(2,1)^2 + v_crf(3,1)^2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% units vector of orbital frame 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% radial 
er = (1 / r_length) * r_crf; 
 
% along-track or tangential component 
et = (1 / v_length) * v_crf; 
 
% cross product 
[cp] = productcross(r_crf,v_crf); 
cp_length = sqrt(cp(1,1)^2 + cp(2,1)^2 + cp(3,1)^2); 
 
% cross-track or normal component 
en = (1 / cp_length) * cp; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Tranformation matrix from inertial to orbital frame 
Rrtn = [er' 
        et' 
        en']; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
     
