function [Xmatrix, Amatrix_rangerate, Wmatrix_rangerate, Amatrix_range, Wmatrix_range, NEQ_range, NEQ_rangerate] = estimator_orbit_intersat(orb1,orb2, veqZarray1,veqParray1, veqZarray2,veqParray2, intersat_obs, sigma_obs) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit parameter estimation based on pseudo-observations and
% intersatellite ranging observations applied as range or range-rate
% constraints
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
% - range_constraint: Baseline length (m) between the 2 satellites applied 
%                     as range constraint in the least squares orbit
%                     estimation of the leading satellite 2
%
% Output arguments:
% - Xmatrix:        Estimated initial sate vector vector of orbit2
% - Amatrix:        Desing matrix per epoch
% - Wmatrix:        b matrix of least squares method  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas D. Papanikolaou                                8 November 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 14/10/2022    Thomas Loudis Papanikolaou 
%               Modified for applying the algorithm to real GRACE-FO intersatellite ranging data 
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
%[dstn,rms_stn,dorb,rms_orb,delta_kepler,rms_kepler,delta_Vstn rms_Vstn] = orbital_pert(orb1,orb2,0);
[dorb,rms_orb,sr] = compstat(orb1,orb2);
%[sz1 sz2] = size(dorb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common epochs between orbits (orb1,orb2) and intersatellite observations (range-rate,range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test time argument of orb2 and intersat_obs
[d1 d2] = size(orbref);
[d3 d4] = size(intersat_obs);
[d5 d6] = size(veqParray2);

Nepochs_orbit = d1;
Nepochs_obs   = d3;
Nparam_veqP   = d6 - Ntime_col;

% Matrices preallocation
% dr_matrix        = zeros(Nepochs,2);
% rangerate_matrix = zeros(Nepochs,2);
% range_matrix     = zeros(Nepochs,3);
% dRTN_matrix      = zeros(Nepochs,4);

iobs0   = 1;
Nepochs = 0;
i_epochs_common = 0;
for iref = 1 : Nepochs_orbit
    % Reference Orbit
    tiref = orbref(iref,1);
    
    % Intersatellite ranging Observations (LRI, KBR data)
    for iobs = iobs0 : Nepochs_obs        
        tiobs = intersat_obs(iobs,1);
        
        % Common Epochs
        if abs(tiref - tiobs) < 10^-8
            i_epochs_common = i_epochs_common + 1;
            ObsEpochs(i_epochs_common,1) = orbref(iref,1);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y observations (measurements)  at common epoch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y measurements :: Range 
y_meas_range = intersat_obs(iobs, Ntime_col + 1);

% y measurements :: Range-Rate
y_meas_rangerate = intersat_obs(iobs, Ntime_col + 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Range
f_obs_range = function_range;
% Range-Rate
f_obs_rangerate = function_rangerate;
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
% b matrix or W matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range
Wmatrix_range(i_epochs_common, 1) = y_meas_range - f_obs_range;

% Range-rate
Wmatrix_rangerate(i_epochs_common, 1) = y_meas_rangerate - f_obs_rangerate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design Matrix : Amatrix_range and Amatrix_rangerate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% veqZarray for epoch ti (State Vector)
VEQ_Z_ti = veqZarray( (iref-1)*6 + 1 : (iref-1)*6 + Nobsset , (Ntime_col + 1) : (Ntime_col + 6) );

% veqParray for epoch ti (Force parameters)
VEQ_P_ti = veqParray( (iref-1)*6 + 1 : (iref-1)*6 + Nobsset , (Ntime_col + 1) : (Ntime_col + Nparam_veqP) );

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
Amatrix_range_ti = [pdv_range_IC_Zo pdv_range_IC_P];
Amatrix_range(i_epochs_common,:) = Amatrix_range_ti;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range-rate
% d_rdot/dXo d_rdot/dYo d_rdot/dZo d_rdot/d_Uxo d_rdot/d_Uyo d_rdot/d_Uzo

% Initial State Vector
pdv_rangerate_IC_Zo = pdv_rangerate_Z_t * VEQ_Z_ti;
% Force Parameters
pdv_rangerate_IC_P  = pdv_rangerate_Z_t * VEQ_P_ti;

% Design matrix for epoch ti 
Amatrix_rangerate_ti = [pdv_rangerate_IC_Zo pdv_rangerate_IC_P];
Amatrix_rangerate(i_epochs_common,:) = Amatrix_rangerate_ti;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD Code
if 1 < 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design Matrix :: range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % values at current epoch (tiref) 
     Amatrix_ti_Xo = pdv_dr_xB * veqZmatrix2(1,1) + ... 
                     pdv_dr_yB * veqZmatrix2(2,1) + ...
                     pdv_dr_zB * veqZmatrix2(3,1); 
     Amatrix_ti_Yo = pdv_dr_xB * veqZmatrix2(1,2) + ... 
                     pdv_dr_yB * veqZmatrix2(2,2) + ... 
                     pdv_dr_zB * veqZmatrix2(3,2);  
     Amatrix_ti_Zo = pdv_dr_xB * veqZmatrix2(1,3) + ... 
                     pdv_dr_yB * veqZmatrix2(2,3) + ... 
                     pdv_dr_zB * veqZmatrix2(3,3); 
     Amatrix_ti_Vxo= pdv_dr_xB * veqZmatrix2(1,4) + ... 
                     pdv_dr_yB * veqZmatrix2(2,4) + ... 
                     pdv_dr_zB * veqZmatrix2(3,4); 
     Amatrix_ti_Vyo= pdv_dr_xB * veqZmatrix2(1,5) + ... 
                     pdv_dr_yB * veqZmatrix2(2,5) + ... 
                     pdv_dr_zB * veqZmatrix2(3,5); 
     Amatrix_ti_Vzo= pdv_dr_xB * veqZmatrix2(1,6) + ... 
                     pdv_dr_yB * veqZmatrix2(2,6) + ... 
                     pdv_dr_zB * veqZmatrix2(3,6); 

     Amatrix_ti =   [Amatrix_ti_Xo Amatrix_ti_Yo Amatrix_ti_Zo Amatrix_ti_Vxo Amatrix_ti_Vyo Amatrix_ti_Vzo];                                   
     Amatrix(iepoch,:) = Amatrix_ti;
     Amatrix_time(iepoch,:) = [tiref Amatrix_ti];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Normal Equations (NEQ)
     %Amatrix_ti_transpose = Amatrix_ti'
     %Amatrix_multiply = Amatrix_ti' * Amatrix_ti

     %NEQn = Amatrix_ti' * Amatrix_ti + NEQn;
     %NEQu = Amatrix_ti' * Wmatrix_ti + NEQu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Design Matrix :: range-rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % values at current epoch (tiref) 
     Amatrix_ti_Xo = pdv_rangerate_xB * veqZmatrix2(1,1) + ... 
                     pdv_rangerate_yB * veqZmatrix2(2,1) + ...
                     pdv_rangerate_zB * veqZmatrix2(3,1) + ...
                     pdv_rangerate_VxB * veqZmatrix2(4,1) + ... 
                     pdv_rangerate_VyB * veqZmatrix2(5,1) + ... 
                     pdv_rangerate_VzB * veqZmatrix2(6,1) ; 

     Amatrix_ti_Yo = pdv_rangerate_xB * veqZmatrix2(1,2) + ... 
                     pdv_rangerate_yB * veqZmatrix2(2,2) + ...
                     pdv_rangerate_zB * veqZmatrix2(3,2) + ...
                     pdv_rangerate_VxB * veqZmatrix2(4,2) + ... 
                     pdv_rangerate_VyB * veqZmatrix2(5,2) + ... 
                     pdv_rangerate_VzB * veqZmatrix2(6,2) ; 

     Amatrix_ti_Zo = pdv_rangerate_xB * veqZmatrix2(1,3) + ... 
                     pdv_rangerate_yB * veqZmatrix2(2,3) + ...
                     pdv_rangerate_zB * veqZmatrix2(3,3) + ...
                     pdv_rangerate_VxB * veqZmatrix2(4,3) + ... 
                     pdv_rangerate_VyB * veqZmatrix2(5,3) + ... 
                     pdv_rangerate_VzB * veqZmatrix2(6,3) ; 
                 
     Amatrix_ti_Vxo= pdv_rangerate_xB * veqZmatrix2(1,4) + ... 
                     pdv_rangerate_yB * veqZmatrix2(2,4) + ...
                     pdv_rangerate_zB * veqZmatrix2(3,4) + ...
                     pdv_rangerate_VxB * veqZmatrix2(4,4) + ... 
                     pdv_rangerate_VyB * veqZmatrix2(5,4) + ... 
                     pdv_rangerate_VzB * veqZmatrix2(6,4) ; 

     Amatrix_ti_Vyo= pdv_rangerate_xB * veqZmatrix2(1,5) + ... 
                     pdv_rangerate_yB * veqZmatrix2(2,5) + ...
                     pdv_rangerate_zB * veqZmatrix2(3,5) + ...
                     pdv_rangerate_VxB * veqZmatrix2(4,5) + ... 
                     pdv_rangerate_VyB * veqZmatrix2(5,5) + ... 
                     pdv_rangerate_VzB * veqZmatrix2(6,5) ; 
                 
     Amatrix_ti_Vzo= pdv_rangerate_xB * veqZmatrix2(1,6) + ... 
                     pdv_rangerate_yB * veqZmatrix2(2,6) + ...
                     pdv_rangerate_zB * veqZmatrix2(3,6) + ...
                     pdv_rangerate_VxB * veqZmatrix2(4,6) + ... 
                     pdv_rangerate_VyB * veqZmatrix2(5,6) + ... 
                     pdv_rangerate_VzB * veqZmatrix2(6,6) ;                  
                                 
     Amatrix_ti =   [Amatrix_ti_Xo Amatrix_ti_Yo Amatrix_ti_Zo Amatrix_ti_Vxo Amatrix_ti_Vyo Amatrix_ti_Vzo];                                   
     Amatrix_rangerate(iepoch,:) = Amatrix_ti;
     %Amatrix_time_rangerate(iepoch,:) = [tiref Amatrix_ti];
end
% End of OLD Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            iobs0 = iobs + 1;
            clear tiref tiobs            
            break
        end
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least-Squares solution & NEQ matrices :: Range
sigma_obs_range = sigma_obs(:,1);
[Xmatrix_range, NEQn_range, NEQu_range] = estimator_neq_sol(Amatrix_range, Wmatrix_range, sigma_obs_range);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least-Squares solution & NEQ matrices :: Range rate
sigma_obs_rangerate = sigma_obs(:,2);
[Xmatrix_rangerate, NEQn_rangerate, NEQu_rangerate] = estimator_neq_sol(Amatrix_rangerate, Wmatrix_rangerate, sigma_obs_rangerate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dXmatrix = Xmatrix_rangerate - Xmatrix_range;
%delta_X_range_rangerate = dXmatrix(1:6,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overall matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEQ matrices
NEQ_range(:,:,1) = NEQn_range;
NEQ_range(:,1,2) = NEQu_range;

NEQ_rangerate(:,:,1) = NEQn_rangerate;
NEQ_rangerate(:,1,2) = NEQu_rangerate;

% A, b matrices
Amatrix = [Amatrix_range; Amatrix_rangerate];
bmatrix = [Wmatrix_range; Wmatrix_rangerate];

% Cv matrix
[n_sigma d2] = size(sigma_obs);
if n_sigma > 1
Cv_obs_sum = [sigma_obs_range; sigma_obs_rangerate];
else
Cv_obs_sum = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least-Squares solution
[Xmatrix, NEQn, NEQu] = estimator_neq_sol(Amatrix, bmatrix, Cv_obs_sum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
