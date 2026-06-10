function [Amatrix_rangerate_out, Amatrix_range_out] = neq_intersat_Amatrix(orb1,orb2, veqZarray1,veqParray1, veqZarray2,veqParray2, intersat_ranging_epochs)  

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
% Normal Equation matrices formulation for intersatellite ranging observations  
% applied for range and range-rate observations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - orb1:           Orbit of the trailing satelllite  
%                   orb1_i= [t r' v']  
%                   t: MJD incdluding fraction of the day  
%                   r,v: Position and Velcotiy vector cartesian coordinates  
% - orb2:           Orbit of the leading satelllite  
%                   orb1_i= [t r' v']  
%                   t: MJD incdluding fraction of the day  
%                   r,v: Position and Velcotiy vector cartesian coordinates  
% - veqZarray1:     Variational Equations solution, State transition matrix 
%                   per epoch for satellite orbit 1 
% - veqZarray2:     Variational Equations solution, State transition matrix 
%                   per epoch for satellite orbit 2 
% - intersat_ranging_epochs:   Inter-Satellite Ranging Epochs  
% 
% Output arguments: 
% - Amatrix:        Desing matrix per epoch 
% - Wmatrix:        b matrix of least squares method   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas D. Papanikolaou                                8 November 2019 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Last modified: 
% 14/10/2022    Thomas Loudis Papanikolaou  
%               Modified for applying the algorithm to real GRACE-FO intersatellite ranging data  
% 11/08/2025    TLP, Modified and renamed to neq_intersat_ranging  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Number of time arguments 
Ntime_col = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% obstype  :  1/2  [r/v] 
obstype = 2; 
if obstype == 1 
    Nobsset = 3; 
elseif obstype == 2 
    Nobsset = 6; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Reference Orbit 
orbref = orb1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Orbit differences: XYZ 
% [dstn,rms_stn,dorb,rms_orb,delta_kepler,rms_kepler,delta_Vstn rms_Vstn] = orbital_pert(orb1,orb2,0); 
% [dorb,rms_orb,sr] = compstat(orb1,orb2); 
% [sz1 sz2] = size(dorb); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Common Epochs of the Observations/Meausurements  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% [delta_series,rms_dsr,series_common] = compstat(orbref,orbobs); 
[epochs_common] = dataseries_commonepochs(orbref,intersat_ranging_epochs); 
% [epochs_common, sr1_common, sr2_common] = dataseries_commonepochs(orbref,intersat_ranging_epochs); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Common epochs between orbits (orb1,orb2) and intersatellite observations (range-rate,range) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Test time argument of orb2 and intersat_ranging_epochs 
[d1 d2] = size(orbref); 
[d3 d4] = size(intersat_ranging_epochs); 
[d5 d6] = size(veqParray2); 
 
Nepochs_orbit = d1; 
Nepochs_obs   = d3; 
Nparam_veqP   = d6 - Ntime_col; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Matrices preallocation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% dr_matrix        = zeros(Nepochs,2); 
% rangerate_matrix = zeros(Nepochs,2); 
% range_matrix     = zeros(Nepochs,3); 
% dRTN_matrix      = zeros(Nepochs,4); 
 
if size(veqParray2,1) > 1  
Nparam = Nparam_veqP; 
else 
Nparam = 0; 
end 
Amatrix_range       = zeros(Nepochs_orbit,6+Nparam); 
Amatrix_rangerate   = zeros(Nepochs_orbit,6+Nparam); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
iobs0   = 1; 
Nepochs = 0; 
i_epochs_common = 0; 
for iref = 1 : Nepochs_orbit 
    % Reference Orbit 
    tiref = orbref(iref,1); 
     
    % Intersatellite ranging Observations (LRI, KBR data) 
    for iobs = iobs0 : Nepochs_obs         
        tiobs = intersat_ranging_epochs(iobs,1); 
         
        % Common Epochs 
        if abs(tiref - tiobs) < 10^-8 
            i_epochs_common = i_epochs_common + 1; 
            % ObsEpochs(i_epochs_common,1) = orbref(iref,1); 
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Observation function values at common epoch 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% State Vector coordinates 
rB(1,1) = orb2(iref,2); 
rB(2,1) = orb2(iref,3); 
rB(3,1) = orb2(iref,4); 
vB(1,1) = orb2(iref,5); 
vB(2,1) = orb2(iref,6); 
vB(3,1) = orb2(iref,7); 
 
rA(1,1) = orb1(iref,2); 
rA(2,1) = orb1(iref,3); 
rA(3,1) = orb1(iref,4); 
vA(1,1) = orb1(iref,5); 
vA(2,1) = orb1(iref,6); 
vA(3,1) = orb1(iref,7); 
     
% rA,rB 
xA = orb1(iref,2); 
yA = orb1(iref,3); 
zA = orb1(iref,4); 
xB = orb2(iref,2); 
yB = orb2(iref,3); 
zB = orb2(iref,4); 
% vA,vB 
VxA = orb1(iref,5); 
VyA = orb1(iref,6); 
VzA = orb1(iref,7); 
VxB = orb2(iref,5); 
VyB = orb2(iref,6); 
VzB = orb2(iref,7); 
 
delta_r_orb12 = rB - rA; 
delta_v_orb12 = vB - vA; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Range :: Satellites distance via orbit differences XYZ 
function_range = sqrt( delta_r_orb12(1,1)^2 + delta_r_orb12(2,1)^2 + delta_r_orb12(3,1)^2 ) ; 
 
rab_vec  = delta_r_orb12 ; 
rab_magn = function_range; 
 
% Line-Of-Sight Vector 
eab_vec = (1 / rab_magn) * rab_vec; 
vab_vec = delta_v_orb12; 
 
% Range-Rate via orbits :: vab_vec' * eab_vec; 
function_rangerate = (vab_vec' * eab_vec); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Partial Derivatives 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Partial derivatives of orbit 2 (to be estimated) w.r.t. initial state vector :: Variational Equations state transition matrix 
veqZarray = veqZarray2; 
veqParray = veqParray2; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% range at ti 
    dr_ti = function_range; 
% range-rate at ti 
    rangerate_ti = function_rangerate; 
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Partial derivatives of range w.r.t. position vector of orbit 2 at current epoch  
    pdv_dr_xB = (xB - xA) / dr_ti; 
    pdv_dr_yB = (yB - yA) / dr_ti; 
    pdv_dr_zB = (zB - zA) / dr_ti; 
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Partial derivatives of range-rate w.r.t. state vector of orbit 2 at current epoch 
    %dr_AB = vec_rB - vec_rA; 
    %dv_AB = vec_vB - vec_vA; 
    dr_AB = rB - rA; 
    dv_AB = vB - vA; 
    dot_dvdr = dv_AB' * dr_AB; 
    % w.r.t. position vector 
    pdv_rangerate_xB = (-1/dr_ti^2) * pdv_dr_xB * dot_dvdr + (VxB-VxA)/dr_ti; 
    pdv_rangerate_yB = (-1/dr_ti^2) * pdv_dr_yB * dot_dvdr + (VyB-VyA)/dr_ti; 
    pdv_rangerate_zB = (-1/dr_ti^2) * pdv_dr_zB * dot_dvdr + (VzB-VzA)/dr_ti; 
    % w.r.t. velocity vector 
    pdv_rangerate_VxB = pdv_dr_xB; 
    pdv_rangerate_VyB = pdv_dr_yB; 
    pdv_rangerate_VzB = pdv_dr_zB;     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% Partials of Range observations w.r.t. State Vector at epoch t (matrix 1x6): 
pdv_range_Z_t = [pdv_dr_xB pdv_dr_yB pdv_dr_zB 0 0 0];  
 
% Partials of Range-Rate observations w.r.t. State Vector at epoch t (matrix 1x6): 
pdv_rangerate_Z_t = [pdv_rangerate_xB pdv_rangerate_yB pdv_rangerate_zB pdv_rangerate_VxB pdv_rangerate_VyB pdv_rangerate_VzB];  
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Design Matrix : Amatrix_range and Amatrix_rangerate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% veqZarray for epoch ti (State Vector) 
VEQ_Z_ti = veqZarray( (iref-1)*6 + 1 : (iref-1)*6 + Nobsset , (Ntime_col + 1) : (Ntime_col + 6) ); 
 
% veqParray for epoch ti (Force parameters) 
if size(veqParray,1) > 1  
VEQ_P_ti = veqParray( (iref-1)*6 + 1 : (iref-1)*6 + Nobsset , (Ntime_col + 1) : (Ntime_col + Nparam_veqP) ); 
else 
VEQ_P_ti = 0; 
end 
 
% Partials of intersat observations w.r.t. state vector at ti (Zti) 
 
% Partials of Range-Rate w.r.t. State Vector at epoch t (matrix 1x6): 
%pdv_rangerate_Z_t = [pdv_rangerate_xB pdv_rangerate_yB pdv_rangerate_zB pdv_rangerate_VxB pdv_rangerate_VyB pdv_rangerate_VzB];  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Range 
% Initial State Vector 
pdv_range_IC_Zo = pdv_range_Z_t * VEQ_Z_ti; 
% Force Parameters 
pdv_range_IC_P  = pdv_range_Z_t * VEQ_P_ti; 
 
% Design matrix for epoch ti  
if size(veqParray,1) > 1  
Amatrix_range_ti = [pdv_range_IC_Zo pdv_range_IC_P]; 
else 
Amatrix_range_ti = pdv_range_IC_Zo; 
end     
Amatrix_range(i_epochs_common,:) = Amatrix_range_ti; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Range-rate 
% d_rdot/dXo d_rdot/dYo d_rdot/dZo d_rdot/d_Uxo d_rdot/d_Uyo d_rdot/d_Uzo 
 
% Initial State Vector 
pdv_rangerate_IC_Zo = pdv_rangerate_Z_t * VEQ_Z_ti; 
% Force Parameters 
pdv_rangerate_IC_P  = pdv_rangerate_Z_t * VEQ_P_ti; 
 
% Design matrix for epoch ti  
if size(veqParray,1) > 1  
Amatrix_rangerate_ti = [pdv_rangerate_IC_Zo pdv_rangerate_IC_P]; 
else 
Amatrix_rangerate_ti = pdv_rangerate_IC_Zo; 
end 
Amatrix_rangerate(i_epochs_common,:) = Amatrix_rangerate_ti; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            iobs0 = iobs + 1; 
            break 
        end 
    end     
end 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Clear empty values 
N_epochs_common = i_epochs_common;  
% Amatrix_range       = zeros(Nepochs_orbit,6+Nparam); 
% Amatrix_rangerate   = zeros(Nepochs_orbit,6+Nparam); 
% Wmatrix_range       = zeros(Nepochs_orbit,1); 
% Wmatrix_rangerate   = zeros(Nepochs_orbit,1);   
 
Amatrix_range_out       = Amatrix_range(1:N_epochs_common,:); 
Amatrix_rangerate_out   = Amatrix_rangerate(1:N_epochs_common,:); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
 
 
