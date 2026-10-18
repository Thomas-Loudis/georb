function [Amatrix_out] = grace_comb_designshift(Amatrix_in, Amatrix_param_index, Nparam_obs_rangerate) 

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
% Function: grace_comb_designshift 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - Amatrix_in          : Design matrix  
% 
% Output arguments: 
% - Amatrix_out         : Design matrix with accelerometers scale partials
%                         collumns shifted 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                          18 september 2026 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%-------------------------------------------------------------------------- 
% Parameter partials collumns indexes
%-------------------------------------------------------------------------- 
col_orbit1_1 = Amatrix_param_index.orbit1(1,1);
col_orbit1_2 = Amatrix_param_index.orbit1(1,2);

col_orbit2_1 = Amatrix_param_index.orbit2(1,1);
col_orbit2_2 = Amatrix_param_index.orbit2(1,2);

% Gravity parameters
col_common_1 = Amatrix_param_index.global(1,1);
col_common_2 = Amatrix_param_index.global(1,2);

col_obs_rangerate_1 = Amatrix_param_index.rangerate_emp(1,1);
col_obs_rangerate_2 = Amatrix_param_index.rangerate_emp(1,2);

col_orbit1_ACC_Scale_start  = Amatrix_param_index.grace1_acc_scale(1,1);
col_orbit1_ACC_Scale_end    = Amatrix_param_index.grace1_acc_scale(1,2);

col_orbit2_ACC_Scale_start  = Amatrix_param_index.grace2_acc_scale(1,1);
col_orbit2_ACC_Scale_end    = Amatrix_param_index.grace2_acc_scale(1,2);

col_RangeRateEMP_1 = col_obs_rangerate_1 ;
col_RangeRateEMP_2 = col_obs_rangerate_2 ;

col_GRAV_indx1 = col_common_1 ;
col_GRAV_indx2 = col_common_2 ;
%-------------------------------------------------------------------------- 

% -------------------------------------------------------------------------
% Accelerometers Scale matrix as Common (global) Parameters :: partials collumns shift
% -------------------------------------------------------------------------    
if Nparam_obs_rangerate > 0

A_obs = Amatrix_in;
% [n1,n2] = size(A_obs)
Amatrix = [ A_obs(:,1 : col_orbit1_ACC_Scale_start-1) A_obs(:, col_orbit2_1 : col_orbit2_ACC_Scale_start-1 ) ... 
            A_obs(:, col_RangeRateEMP_1 : col_RangeRateEMP_2 ) ...  
            A_obs(:, col_orbit1_ACC_Scale_start : col_orbit1_ACC_Scale_end )  A_obs(:, col_orbit2_ACC_Scale_start : col_orbit2_ACC_Scale_end ) ...
            A_obs(:, col_GRAV_indx1 : col_GRAV_indx2 ) ];
Amatrix_out = Amatrix;

else

A_obs = Amatrix_in;
% [n1,n2] = size(A_obs)
Amatrix = [ A_obs(:,1 : col_orbit1_ACC_Scale_start-1) A_obs(:, col_orbit2_1 : col_orbit2_ACC_Scale_start-1 ) ... 
            A_obs(:, col_orbit1_ACC_Scale_start : col_orbit1_ACC_Scale_end )  A_obs(:, col_orbit2_ACC_Scale_start : col_orbit2_ACC_Scale_end ) ...
            A_obs(:, col_GRAV_indx1 : col_GRAV_indx2 ) ];
Amatrix_out = Amatrix;

end

