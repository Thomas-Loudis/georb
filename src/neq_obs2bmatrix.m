function [bmatrix] = neq_obs2bmatrix (orbref,orbobs) 

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
% b matrix: Reduced Observation matrix (Measurements minus Observation function values) 
% Formulated at the common epochs between the pseudo-Observation matrix and 
% the computed dynamic orbit  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - orbref:         Dynamic Orbit (observation function) at the common 
%                   epochs with the observations (measurements) 
%                   orbref = [t r' v']  
%                   t: MJD incdluding fraction of the day  
%                   r,v: Position and Velcotiy vector cartesian coordinates  
% - orbobs:         Pseudo-Observations based on Kinematic Orbit positions (XYZ) 
%                   at the common epochs with the dynamic orbit  
%                   orbobs = [t r'] or orbobs = [t r' v']  
%                   t: MJD incdluding fraction of the day  
%                   r,v: Position and Velcotiy vector cartesian coordinates  
% 
% Output arguments: 
% - bmatrix:        b matrix Reduced Observation matrix (Measurements minus Observation function values)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                             11 August 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
% [n1_orbref,n2_orbref] = size(orbref); 
% [n1_orbobs,n2_orbobs] = size(orbobs); 
 
orbref_common = orbref;  
orbobs_common = orbobs;  
epochs_common = orbref(:,1); 
 
% obstype  :  1/2  [r/v] 
obstype = 1; 
if obstype == 1 
    Nobsset = 3; 
elseif obstype == 2 
    Nobsset = 6; 
end 
 
% Matrices Preallocation 
[N_commonepochs, n2] = size(epochs_common); 
Wmatrix = zeros(N_commonepochs * Nobsset, 1); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% b matrix :: Reduced Observation matrix (Meaurements minus Observation functions values) formed for the common epochs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for i_epochs = 1 : N_commonepochs 
    iref = i_epochs; 
    iobs = i_epochs; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    Wmatrix_ti = [ orbobs_common(iobs,2:(Nobsset+1)) - orbref_common(iref,2:(Nobsset+1)) ]'; 
    Wmatrix((i_epochs-1)*Nobsset+1 : i_epochs*Nobsset,1) = Wmatrix_ti; 
    % Wmatrix_time((i_epochs-1)*Nobsset+1 : i_epochs*Nobsset,:) = [timearray Wmatrix_ti]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end 
 
bmatrix = Wmatrix; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
