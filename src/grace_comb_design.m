function [A_orbit1, A_orbit2, A_range, A_rangerate, Amatrix_param_index] = grace_comb_design(design_matrix, N_param_GRAV, Nglobal_scale, Nparam_rangerate) 

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
% Function: grace_comb_design 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - design_matrix       : Design matrix structure array including all Amatrices 
% 
% Output arguments: 
% - design_matrix       : Design matrix structure array with processed design matrices 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                          18 september 2026 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%-------------------------------------------------------------------------- 
% Design matrix matrices
%-------------------------------------------------------------------------- 
Amatrix_rangerate_sat1 = design_matrix.rangerate_sat1;
Amatrix_range_sat1     = design_matrix.range_sat1;
Amatrix_rangerate_sat2 = design_matrix.rangerate_sat2;
Amatrix_range_sat2     = design_matrix.range_sat2;
Amatrix_obsorb_grace1  = design_matrix.obsorb_grace1;
Amatrix_obsorb_grace2  = design_matrix.obsorb_grace2;
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Parameter types numbers
%-------------------------------------------------------------------------- 
Nparam_common = N_param_GRAV;
N_common_ACC_Scale = Nglobal_scale;

% Range-Rate Empirical parameters 
Nparam_obs_rangerate = Nparam_rangerate;
%-------------------------------------------------------------------------- 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dimensions of design matrix and Number of parameters, observations  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[d1, d2] = size(Amatrix_obsorb_grace1); 
[d3, d4] = size(Amatrix_obsorb_grace2); 
Nparam_grace1 = d2; 
Nparam_grace2 = d4; 
Nobs_sat1 = d1; 
Nobs_sat2 = d3;

% Overall parameters number of the combined parameter estimation including
% gravity parameters (global parameters) 
% Nparam_COMB = (Nparam_grace1 - N_param_GRAV) + (Nparam_grace2 - N_param_GRAV) + N_param_GRAV + Nparam_obs_rangerate;
Nparam_COMB = Nparam_grace1 + Nparam_grace2 - Nparam_common + Nparam_obs_rangerate; 
% Number of orbital parameters (individual arc parameters) 
Nparam_orbital = Nparam_COMB - Nparam_common; 
 
[d5, d6] = size(Amatrix_rangerate_sat1); 
[d7, d8] = size(Amatrix_rangerate_sat2); 
Nobs_rangerate = d7; 
Nparam_rangerate = d8; 
 
[ d9, d10] = size(Amatrix_range_sat1); 
[d11, d12] = size(Amatrix_range_sat2); 
Nobs_range = d9; 
Nparam_range = d10; 
 
N_obs_comb = Nobs_sat1 + Nobs_sat2 + Nobs_rangerate + Nobs_range; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Design Matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Design matrices' Collumns index numbers 
col_orbit1_1 = 1 ;
col_orbit1_2 = Nparam_grace1-Nparam_common ;
 
col_orbit2_1 = col_orbit1_2 + 1 ; 
col_orbit2_2 = col_orbit1_2 + Nparam_grace2-Nparam_common ; 
 
col_obs_rangerate_1 = col_orbit2_2 + Nparam_obs_rangerate / Nparam_obs_rangerate ; 
col_obs_rangerate_2 = col_orbit2_2 + Nparam_obs_rangerate; 
 
col_common_1 = col_orbit2_2 + Nparam_obs_rangerate + 1 ;
col_common_2 = col_orbit2_2 + Nparam_obs_rangerate + Nparam_common ;
 
% Selected parameters to be estimated :: All or Global parameters (common) 
param_col_start = 1; 
param_col_end = Nparam_COMB; 

%--------------------------------------------------------------------------  
col_RangeRateEMP_1 = col_obs_rangerate_1 ;
col_RangeRateEMP_2 = col_obs_rangerate_2 ;

col_orbit1_ACC_Scale_start = col_orbit1_2 - N_common_ACC_Scale/2 + 1 ;
col_orbit1_ACC_Scale_end   = col_orbit1_2 ;

col_orbit2_ACC_Scale_start = col_orbit2_2 - N_common_ACC_Scale/2 + 1 ;
col_orbit2_ACC_Scale_end   = col_orbit2_2 ;

col_GRAV_indx1 = Nparam_orbital + 1 ;
col_GRAV_indx2 = Nparam_orbital + N_param_GRAV ;
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Design matrix :: Parameters collumn indexes
Amatrix_param_index.orbit1 = [col_orbit1_1 col_orbit1_2];
Amatrix_param_index.orbit2 = [col_orbit2_1 col_orbit2_2];
Amatrix_param_index.global = [col_common_1 col_common_2];
Amatrix_param_index.rangerate_emp    = [col_obs_rangerate_1 col_obs_rangerate_2];
Amatrix_param_index.grace1_acc_scale = [col_orbit1_ACC_Scale_start col_orbit1_ACC_Scale_end];
Amatrix_param_index.grace2_acc_scale = [col_orbit2_ACC_Scale_start col_orbit2_ACC_Scale_end];
%--------------------------------------------------------------------------


%-------------------------------------------------------------------------- 
% Orbit 1 
A_orbit1 = zeros(Nobs_sat1 , param_col_end - (param_col_start-1)); 
% Orbital Elements 
A_orbit1(1 : Nobs_sat1 , col_orbit1_1 : col_orbit1_2) = Amatrix_obsorb_grace1(:,1: col_orbit1_2); 
% Common parameters :: Gravity 
if Nparam_common > 0 
A_orbit1(1 : Nobs_sat1 , col_common_1 : col_common_2) = Amatrix_obsorb_grace1(:, Nparam_grace1-Nparam_common+1 : Nparam_grace1-Nparam_common+Nparam_common); 
end 
 
% Orbit 2 
A_orbit2 = zeros(Nobs_sat2 , param_col_end - (param_col_start-1)); 
A_orbit2(1 : Nobs_sat2 , col_orbit2_1 : col_orbit2_2) = Amatrix_obsorb_grace2(:,1:Nparam_grace2-Nparam_common); 
if Nparam_common > 0 
A_orbit2(1 : Nobs_sat2 , col_common_1 : col_common_2) = Amatrix_obsorb_grace2(:,Nparam_grace2-Nparam_common+1 : Nparam_grace2-Nparam_common+Nparam_common); 
end 

intersat_ranging_partials_sat12 = 1;

% Range-rate 
col_common_0 = Nparam_rangerate - Nparam_common - Nparam_obs_rangerate;
if intersat_ranging_partials_sat12 == 1
A_rangerate(1 : Nobs_rangerate , col_orbit1_1 : col_orbit1_2) = Amatrix_rangerate_sat1(:,1:col_common_0); 
end
A_rangerate(1 : Nobs_rangerate , col_orbit2_1 : col_orbit2_2) = Amatrix_rangerate_sat2(:,1:col_common_0); 
% Observation Empirical Parameters 
if Nparam_obs_rangerate > 0 
% A_rangerate(1 : Nobs_rangerate , col_obs_rangerate_1 : col_obs_rangerate_2) = Amatrix_rangerate_sat1(:,col_common_0+Nparam_common+1:col_common_0+Nparam_common+Nparam_obs_rangerate); % + Amatrix_rangerate_sat2(:,col_common_0+Nparam_common+1:col_common_0+Nparam_common+Nparam_obs_rangerate);  
A_rangerate(1 : Nobs_rangerate , col_obs_rangerate_1 : col_obs_rangerate_2) = Amatrix_rangerate_sat2(:,col_common_0+Nparam_common+1:col_common_0+Nparam_common+Nparam_obs_rangerate);  
end 
% Common parameters :: Gravity 
if Nparam_common > 0 
A_rangerate(1 : Nobs_rangerate , col_common_1 : col_common_2) = Amatrix_rangerate_sat1(:,col_common_0+1:col_common_0+Nparam_common) + Amatrix_rangerate_sat2(:,col_common_0+1:col_common_0+Nparam_common);  
% A_rangerate(1 : Nobs_rangerate , col_common_1 : col_common_2) = Amatrix_rangerate_sat1(:,col_common_0+1:col_common_0+Nparam_common); % + Amatrix_rangerate_sat2(:,col_common_0+1:col_common_0+Nparam_common);  
% A_rangerate(1 : Nobs_rangerate , col_common_1 : col_common_2) = Amatrix_rangerate_sat2(:,col_common_0+1:col_common_0+Nparam_common);  
end 
 
% Range 
col_common_0 = Nparam_range - Nparam_common - Nparam_obs_rangerate; 
if intersat_ranging_partials_sat12 == 1
A_range(1 : Nobs_range, col_orbit1_1 : col_orbit1_2) = Amatrix_range_sat1(:,1:col_common_0); 
end
A_range(1 : Nobs_range, col_orbit2_1 : col_orbit2_2) = Amatrix_range_sat2(:,1:col_common_0); 
if Nparam_obs_rangerate > 0 
% A_range(1 : Nobs_range, col_obs_rangerate_1 : col_obs_rangerate_2) = Amatrix_range_sat1(:,col_common_0+Nparam_common+1:col_common_0+Nparam_common+Nparam_obs_rangerate); % + Amatrix_range_sat2(:,col_common_0+1:col_common_0+Nparam_obs_rangerate); 
A_range(1 : Nobs_range, col_obs_rangerate_1 : col_obs_rangerate_2) = Amatrix_range_sat2(:,col_common_0+Nparam_common+1:col_common_0+Nparam_common+Nparam_obs_rangerate); 
end 
if Nparam_common > 0 
A_range(1 : Nobs_range, col_common_1 : col_common_2) = Amatrix_range_sat1(:,col_common_0+1:col_common_0+Nparam_common) + Amatrix_range_sat2(:,col_common_0+1:col_common_0+Nparam_common); 
% A_range(1 : Nobs_range, col_common_1 : col_common_2) = Amatrix_range_sat1(:,col_common_0+1:col_common_0+Nparam_common); % + Amatrix_range_sat2(:,col_common_0+1:col_common_0+Nparam_common); 
% A_range(1 : Nobs_range, col_common_1 : col_common_2) = Amatrix_range_sat2(:,col_common_0+1:col_common_0+Nparam_common); 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% if NEQ_write == 1 
% save A_orbit1.neq  A_orbit1 -ASCII -double 
% save A_orbit2.neq  A_orbit2 -ASCII -double 
% save A_rangerate.neq  A_rangerate -ASCII -double 
% save A_range.neq  A_range -ASCII -double 
% 
% save b_orbit1.neq  b_orbit1 -ASCII -double 
% save b_orbit2.neq  b_orbit2 -ASCII -double 
% save b_rangerate.neq  b_rangerate -ASCII -double 
% save b_range.neq  b_range -ASCII -double 
% end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
% -------------------------------------------------------------------------
% Accelerometers Scale matrix as Common (global) Parameters :: partials collumns shift
% -------------------------------------------------------------------------    
if N_common_ACC_Scale > 0
    
[A_orbit1]      = grace_comb_designshift(A_orbit1, Amatrix_param_index, Nparam_obs_rangerate);
[A_orbit2]      = grace_comb_designshift(A_orbit2, Amatrix_param_index, Nparam_obs_rangerate);
[A_rangerate]   = grace_comb_designshift(A_rangerate, Amatrix_param_index, Nparam_obs_rangerate);
[A_range]       = grace_comb_designshift(A_range, Amatrix_param_index, Nparam_obs_rangerate);

%-------------------------------------------------------------------------- 
% Update Parameter Collumn Indexes considering the shift for accelerometers
% parameters as global
%-------------------------------------------------------------------------- 
% col_orbit1_1 = 1 ;
col_orbit1_2 = col_orbit1_2 - N_common_ACC_Scale/2 ;

% col_orbit2_1 = col_orbit2_1 - N_common_ACC_Scale/2 
% col_orbit2_2 = col_orbit2_2 - N_common_ACC_Scale/2 
col_orbit2_1 = col_orbit1_2 + 1 ;
col_orbit2_2 = col_orbit1_2 + Nparam_grace2-Nparam_common - N_common_ACC_Scale/2 ;

% col_obs_rangerate_1 = col_obs_rangerate_1 - N_common_ACC_Scale/2 
% col_obs_rangerate_2 = col_obs_rangerate_2 - N_common_ACC_Scale/2 
col_obs_rangerate_1 = col_orbit2_2 + Nparam_obs_rangerate / Nparam_obs_rangerate ;
col_obs_rangerate_2 = col_orbit2_2 + Nparam_obs_rangerate ; 

col_orbit1_ACC_Scale_start = col_obs_rangerate_2 + 1;
col_orbit1_ACC_Scale_end   = col_obs_rangerate_2 + N_common_ACC_Scale/2;
col_orbit2_ACC_Scale_start = col_orbit1_ACC_Scale_end + 1 ;
col_orbit2_ACC_Scale_end   = col_orbit1_ACC_Scale_end + N_common_ACC_Scale/2 ;
%-------------------------------------------------------------------------- 

end
% -------------------------------------------------------------------------

%-------------------------------------------------------------------------- 
% Combined Design matrices
%-------------------------------------------------------------------------- 
design_matrix.rangerate_comb    = A_rangerate;
design_matrix.range_comb        = A_range;
design_matrix.orbit1_comb       = A_orbit1;
design_matrix.orbit2_comb       = A_orbit2;
%--------------------------------------------------------------------------
% Design matrix :: Parameters collumn indexes
Amatrix_param_index.orbit1 = [col_orbit1_1 col_orbit1_2];
Amatrix_param_index.orbit2 = [col_orbit2_1 col_orbit2_2];
Amatrix_param_index.global = [col_common_1 col_common_2];
Amatrix_param_index.rangerate_emp    = [col_obs_rangerate_1 col_obs_rangerate_2];
Amatrix_param_index.grace1_acc_scale = [col_orbit1_ACC_Scale_start col_orbit1_ACC_Scale_end];
Amatrix_param_index.grace2_acc_scale = [col_orbit2_ACC_Scale_start col_orbit2_ACC_Scale_end];
%--------------------------------------------------------------------------
