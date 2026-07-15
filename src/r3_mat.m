function [R3_theta] = r3_mat(theta) 

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
% Function: r3 
%-------------------------------------------------------------------------- 
% Purpose 
%  Rotation matrix around axis Z 
%-------------------------------------------------------------------------- 
% Input arguments: 
% - theta:         Angle in radians
% 
% Output arguments: 
% - R3_theta:      Rotation matrix around Z axis for theta angle,
%                  considering anticlockwise rotation
%-------------------------------------------------------------------------- 
% Thomas Loudis Papanikolaou                                   15 July 2026 
%-------------------------------------------------------------------------- 
 
 
%-------------------------------------------------------------------------- 
% Rotation matrix R3 or Rz
%-------------------------------------------------------------------------- 
R3_theta = [ cos(theta)   sin(theta)     0 
            -sin(theta)   cos(theta)     0 
                 0             0           1];               
%-------------------------------------------------------------------------- 
