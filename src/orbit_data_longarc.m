function [orbit_model_struct] = orbit_data_longarc (main_config_fname, orbit_model_filename, config_struct, src_version, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_data_longarc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Orbit data over long orbit arc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                  29 April 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cfg_fname = config_struct;

% Orbit model matrix 
IC_MJDo = orbit_model_struct.IC_MJD ;
% Seconds since start of the day (00h) in TT
IC_Sec_00 = (IC_MJDo - fix(IC_MJDo)) * 86400;
% orbit_model_struct.IC_CRF = [IC_MJDo IC_Zo_vec'];
orbit_arc_length_sec = orbit_model_struct.orbit_arc_length_sec;

N_daily_arcs = ceil(orbit_arc_length_sec / 86400);

param_keyword = 'orbiting_object_name';
[orbiting_object_name] = read_param_cfg(cfg_fname,param_keyword);

param_keyword = 'Reference_frame';
[IC_reference_frame] = read_param_cfg(cfg_fname,param_keyword);

% param_keyword = 'Time_scale';
% [IC_time_scale] = read_param_cfg(cfg_fname,param_keyword);

% param_keyword = 'Date_format';
% [IC_epoch_date_format] = read_param_cfg(cfg_fname,param_keyword);

param_keyword = 'N_orbit_arcs';
[N_orbit_arcs_time_series] = read_param_cfg(cfg_fname,param_keyword);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial State Vector
param_keyword = 'State_vector';
[param_value, param_line] = read_param_cfg(cfg_fname,param_keyword);
IC_Zo(1,1) = sscanf(param_line,'%f %*');
IC_Zo(2,1) = sscanf(param_line,'%*s %f %*');
IC_Zo(3,1) = sscanf(param_line,'%*s %*s %f %*');
IC_Zo(4,1) = sscanf(param_line,'%*s %*s %*s %f %*');
IC_Zo(5,1) = sscanf(param_line,'%*s %*s %*s %*s %f %*');
IC_Zo(6,1) = sscanf(param_line,'%*s %*s %*s %*s %*s %f %*');
IC_line = param_line;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_arc = 1 : N_daily_arcs
    IC_time_scale = 'TT';
    IC_epoch_date_format = 'MJD';
    % Update MJD
    MJD_day = fix(IC_MJDo) + (i_arc-1);
    Sec_00  = IC_Sec_00;
    % Write IC data line format
    ic_data_line = sprintf('%s %s %s %s %d %f %s %s %s %s',orbiting_object_name, IC_reference_frame, IC_time_scale, IC_epoch_date_format, MJD_day, Sec_00, '0', '0', N_orbit_arcs_time_series, IC_line); 
    ic_data_array(i_arc,:) = ic_data_line;
end

% IC data array of sequential epochs
[ic_n, ic_m] = size(ic_data_array);
for ic_i = 1 : ic_n
    % IC per ndividual epoch
    ic_data_object_i = ic_data_array(ic_i,:);
    
    % Orbit configuration structure
    [orbit_config_struct] = write_config2struct(main_config_fname, orbit_model_filename, ic_data_object_i, src_version);
    
    % Orbit Data read per IC epoch
    [accelerometer_struct, observation_struct, KBR_intersat_struct, LRI_intersat_struct] = orbit_data (orbit_config_struct, orbit_model_struct);

% Auxiliary matrices for the passing temporarily the individaul arcs data    
if ic_i == 1
acc_data_matrix = accelerometer_struct.ACC1B_data_array;
sca_data_matrix = accelerometer_struct.SCA1B_data_array;

obs_matrix_crf  = observation_struct.obs_orbit_crf;
obs_matrix_trf  = observation_struct.obs_orbit_trf;

kbr_range = KBR_intersat_struct.range;
kbr_rangerate  = KBR_intersat_struct.rangerate;
kbr_rangeaccel = KBR_intersat_struct.rangeacceleration;

lri_range = LRI_intersat_struct.range;
lri_rangerate  = LRI_intersat_struct.rangerate;
lri_rangeaccel = LRI_intersat_struct.rangeacceleration;

else
acc_data_matrix = [acc_data_matrix; accelerometer_struct.ACC1B_data_array];
sca_data_matrix = [sca_data_matrix; accelerometer_struct.SCA1B_data_array];

obs_matrix_crf  = [obs_matrix_crf; observation_struct.obs_orbit_crf];
obs_matrix_trf  = [obs_matrix_trf; observation_struct.obs_orbit_trf];

kbr_range = [kbr_range; KBR_intersat_struct.range] ;
kbr_rangerate  = [kbr_rangerate; KBR_intersat_struct.rangerate] ;
kbr_rangeaccel = [kbr_rangeaccel; KBR_intersat_struct.rangeacceleration] ;

lri_range = [lri_range; LRI_intersat_struct.range] ;
lri_rangerate  = [lri_rangerate; LRI_intersat_struct.rangerate] ;
lri_rangeaccel = [lri_rangeaccel; LRI_intersat_struct.rangeacceleration] ;   
end

end
% End of daily arcs loop

% Updating strucutre arrays with the overall matrices of the long arcs
accelerometer_struct.ACC1B_data_array = acc_data_matrix;
accelerometer_struct.SCA1B_data_array = sca_data_matrix;

observation_struct.obs_orbit_crf = obs_matrix_crf;
observation_struct.obs_orbit_trf = obs_matrix_trf;

KBR_intersat_struct.range = kbr_range;
KBR_intersat_struct.rangerate = kbr_rangerate;
KBR_intersat_struct.rangeacceleration = kbr_rangeaccel;

LRI_intersat_struct.range = lri_range;
LRI_intersat_struct.rangerate = lri_rangerate;
LRI_intersat_struct.rangeacceleration = lri_rangeaccel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update central "orbit model structure matrix"
orbit_model_struct.accelerometer_struct = accelerometer_struct;
orbit_model_struct.observation_struct   = observation_struct;
orbit_model_struct.KBR_intersat_struct  = KBR_intersat_struct;
orbit_model_struct.LRI_intersat_struct  = LRI_intersat_struct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
