function [R2_theta] = r2_mat(theta) 

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
% Function: r2_mat 
%-------------------------------------------------------------------------- 
% Purpose 
%  Rotation matrix around axis X 
%-------------------------------------------------------------------------- 
% Input arguments: 
% - theta:         Angle in radians
% 
% Output arguments: 
% - R2_theta:      Rotation matrix around X axis for theta angle 
%                  considering anticlockwise rotation
%-------------------------------------------------------------------------- 
% Thomas Loudis Papanikolaou                                   15 July 2026 
%-------------------------------------------------------------------------- 
 
 
%-------------------------------------------------------------------------- 
% Rotation matrix R2 or Ry
%-------------------------------------------------------------------------- 
R2_theta = [ cos(theta)    0    -sin(theta) 
                 0         1          0 
             sin(theta)    0     cos(theta) ]; 
%-------------------------------------------------------------------------- 
