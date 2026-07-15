function [R1_theta] = r1_mat(theta) 

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
 
 
%-------------------------------------------------------------------------- 
% Function: r1_mat 
%-------------------------------------------------------------------------- 
% Purpose 
%  Rotation matrix around axis X 
%-------------------------------------------------------------------------- 
% Input arguments: 
% - theta:         Angle in radians
% 
% Output arguments: 
% - R1_theta:      Rotation matrix around X axis for theta angle 
%                  considering anticlockwise rotation
%-------------------------------------------------------------------------- 
% Thomas Loudis Papanikolaou                                   15 July 2026 
%-------------------------------------------------------------------------- 
 
 
%-------------------------------------------------------------------------- 
% Rotation matrix R1 or Rx
%-------------------------------------------------------------------------- 
R1_theta = [   1        0            0 
               0    cos(theta)   sin(theta) 
               0   -sin(theta)   cos(theta) ]; 
%-------------------------------------------------------------------------- 

