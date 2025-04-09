function [orbit_config_struct, orbit_matrix, orbit_rms, veqZarray, veqParray, OBS_matrix, Xparam_aposteriori, OBS_residuals, ext_orb, forces_accel, orbit_model_struct] = orbit_object(orbit_config_struct, write_data, orbit_model_struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  orbit_object performs the individual orbit analysis (precise orbit
%  determination, orbit propagation, simulation) of a single orbiting
%  object (satellite, individual object) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - main_config_fname     : GEORB master configuration file name
% - ic_data_object        : 
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             23 August 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Determination mode
param_keyword = 'orbit_mode';
[orbit_mode] = read_param_cfg(orbit_config_struct,param_keyword);
test = strcmp(orbit_mode,'orbit_propagation_veq');
if test == 1
    fprintf('\n%s\n','Mode: Orbit Propagation EQM and VEQ');    
end
test = strcmp(orbit_mode,'orbit_propagation_eqm');
if test == 1
    fprintf('\n%s\n','Mode: Orbit Propagation EQM');    
end
test = strcmp(orbit_mode,'orbit_determination');
if test == 1
    fprintf('\n%s\n','Mode: Orbit Determination');    
end

% Satellite or Orbiting OBject name
param_keyword = 'orbiting_object_name';
[orbiting_object_name] = read_param_cfg(orbit_config_struct,param_keyword);
fprintf('%s%s\n','Orbiting Object: ', orbiting_object_name);

% Initial Epoch
param_keyword = 'Initial_Epoch';
[param_value, param_line] = read_param_cfg(orbit_config_struct,param_keyword);
Initial_Epoch = param_line;
fprintf('%s%s\n','Initial Epoch: ',Initial_Epoch);
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Determination :: Call POD function orbit_pod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = strcmp(orbit_mode,'orbit_determination');
if test == 1
fprintf('%s\n', 'Orbit Determination :: in-progress : Iterations of Orbit Integration & Parameter Estimation');   
else
fprintf('%s\n', 'Orbit Propagation :: in-progress');   
end   
% fprintf('%s\n','Models and Data preprocessing :: In-progress');    
% Call function orbit_pod
[rms_orbital,rms_orbc,rms_orbk,rms_orbt,dstn,dorbc,dkepl,dorbt,Xmatrix,orbc,orbk,orbt,veqZarray,veqParray,rms_obs,OBS_matrix, Xparam_aposteriori, OBS_residuals, ext_orb, forces_accel, Gmatrix, Rmatrix, orbit_model_struct] = orbit_pod(orbit_config_struct, orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = strcmp(orbit_mode,'orbit_determination');
if test == 1
    % fprintf('\n');
    fprintf('%s %11.6f %11.6f %11.6f', 'Orbit residuals: RMS(XYZ): ',rms_obs(1:3));
    fprintf('\n\n');
end

param_keyword = 'external_orbit_comp';
[external_orbit_comp_yn] = read_param_cfg(orbit_config_struct,param_keyword);
test = strcmp(external_orbit_comp_yn,'y');
if test == 1
    fprintf('%s \n', 'External Orbit Comparison:');
    fprintf('%s ', 'Orbit residuals: RMS(RTN): ' );
    fprintf('%11.6f ', rms_orbital);
    fprintf('\n');
    fprintf('%s ', 'Orbit residuals: RMS(XYZ): ' );
    fprintf('%11.6f ', rms_orbc(1:3));
    fprintf('\n\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbits in 3-dimensional matrix
[Nepochs, Nelements] = size(orbc);
orbit_matrix = zeros(Nepochs, Nelements, 3);
orbit_matrix(:,:,1) = orbc;
orbit_matrix(:,:,2) = orbt;
orbit_matrix(:,:,3) = orbk;

% Statistics matrix of orbit residuals and external orbit comparison
orbit_rms      = zeros(5,6);
orbit_rms(1,1:3) = rms_obs;
orbit_rms(2,1:3) = rms_orbital;
orbit_rms(3,:) = rms_orbc;
orbit_rms(4,:) = rms_orbt;
orbit_rms(5,:) = rms_orbk;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write orbit and partial derivatives to files in georb format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if write_data == 1
%[fid] = write_georb_data(orbit_config_struct, orbc, orbt, orbk, veqZarray, veqParray, OBS_residuals, dstn,dorbc,dkepl);

% Make Directory for moving the data files
[OUT_fname_object_mjd, OUT_fname_mission_mjd] = write_results_dir(orbit_config_struct,orbit_model_struct);
[status, message, messageid] = rmdir(OUT_fname_object_mjd,'s');
[status, message, messageid] = mkdir(OUT_fname_object_mjd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbits
data_matrix = orbc;
data_functional = 'Orbit';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);

data_matrix = orbt;
data_functional = 'Orbit';
reference_frame = 'ITRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);

data_matrix = orbk;
data_functional = 'Orbit';
reference_frame = 'Kepler';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial Derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n1, n2] = size(veqZarray);
[n3, n4] = size(veqParray);
if n1 > 1 
    if n1 == n3
    VEQ_matrix = [veqZarray veqParray(:,2:n4)];
    end
    if n3 <= 1 
    VEQ_matrix = [veqZarray];
    end    
data_matrix = VEQ_matrix;
data_functional = 'Partial Derivatives';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);


if 1 < 0

% VEQ-Z State Transition Matrix
data_matrix = veqZarray;
data_functional = 'Partials_StateTransitionMatrix';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);
% Reshape matrix-based format to row format
No_col_time = 1;
matrix_step = 6;
[data_matrix_2] = matrix_reshape(data_matrix, No_col_time, matrix_step);
% GNP format
[georb2gnv_data_name] = write_georb_data3_gnv(orbit_config_struct, data_matrix_2, data_functional, reference_frame);
[status,message,messageid] = movefile(georb2gnv_data_name, OUT_fname_object_mjd);

if n3 > 1 
% VEQ-P Sensitivity Matrix
data_matrix = veqParray;
data_functional = 'Partials_SensitivityMatrix';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);


% VEQ-P Sensitivity Matrix w.r.t. Accelerometer Calibration Parameters
data_matrix = veqParray(:,1:7);
data_functional = 'Partials_SensitivityMatrix_ACC-CAL';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);
% Reshape matrix-based format to row format
[data_matrix_2] = matrix_reshape(data_matrix, No_col_time, matrix_step);
% GNP format
[georb2gnv_data_name] = write_georb_data3_gnv(orbit_config_struct, data_matrix_2, data_functional, reference_frame);
[status,message,messageid] = movefile(georb2gnv_data_name, OUT_fname_object_mjd);


GRAV_param_estim = 0;

if GRAV_param_estim == 1
% VEQ-P Sensitivity Matrix w.r.t. Gravity Field Parameters
data_matrix(:,1) = veqParray(:,1);
Ngravparam = n4 -1 -6;
orbit_object_writedata_Ngravparam = Ngravparam
% Ntimeargument = 1; Nacc_cal_param = 6;
data_matrix(:, 2 : Ngravparam + 1) = veqParray(:, 8 : Ngravparam + 7);
data_functional = 'Partials_SensitivityMatrix_GRAV';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);
% Reshape matrix-based format to row format
[data_matrix_2] = matrix_reshape(data_matrix, No_col_time, matrix_step);
% GNP format
[georb2gnv_data_name] = write_georb_data3_gnv(orbit_config_struct, data_matrix_2, data_functional, reference_frame);
[status,message,messageid] = movefile(georb2gnv_data_name, OUT_fname_object_mjd);
end

end

end
% end of writting individual partials data files

end
% end of partials 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forces: Acceleration 
data_matrix = forces_accel;
data_functional = 'Forces acceleration';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);

if 1 < 0 
% Forces: Gravity_gradient 
data_matrix = Gmatrix;
data_functional = 'Gravity_Gradient';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);
% GNP format
[georb2gnv_data_name] = write_georb_data3_gnv(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb2gnv_data_name, OUT_fname_object_mjd);

% Earth Rotation: Earth Orientation matrix 
data_matrix = Rmatrix;
data_functional = 'Earth_Rotation';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);
% GNP format
[georb2gnv_data_name] = write_georb_data3_gnv(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb2gnv_data_name, OUT_fname_object_mjd);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit observation residuals
[n1 n2] = size(OBS_residuals);
if n1 > 1 
data_matrix = OBS_residuals;
data_functional = 'Orbit Observations residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters (Vector matrix)
data_matrix = Xparam_aposteriori';
data_functional = 'Parameters';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);

% Parameters in GNP1B data file
[georb2gnv_data_name] = write_georb_data3_gnv(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb2gnv_data_name, OUT_fname_object_mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External Orbit comparison
[n1 n2] = size(dstn);
if n1 > 1 
data_matrix = dstn;
data_functional = 'External Orbit comparison';
reference_frame = 'Orbital Frame';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);

data_matrix = dkepl;
data_functional = 'External Orbit comparison';
reference_frame = 'Kepler';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);

data_matrix = dorbc;
data_functional = 'External Orbit comparison';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Statistics
rms_kbr = 0;
rms_lri = 0;
orbit_config_fname_pair = 'no_pair';
[georb_data_name] = write_georb_statistics(orbit_config_struct, orbit_config_fname_pair, rms_obs, rms_orbital, rms_kbr, rms_lri);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Configuration parameters :: Structure matrix to external config file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output data file name conventions
[georb_dataformat_name, georb_dataformat_suffix] = write_data_name(orbit_config_struct, 0);
% Data functional name
data_functional_filename = '_config';
% Data file name
georb_data_name = sprintf('%s%s%s',georb_dataformat_name,data_functional_filename,georb_dataformat_suffix);
[orbit_config_filename] = write_config_struct2file(orbit_config_struct,georb_data_name);
[status,message,messageid] = movefile(georb_data_name, OUT_fname_object_mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GNP1B Orbits in GNV data format (GRACE missions orbit data format)
mission_name = ' '; 
param_keyword = 'georb_mode';
[georb_mode] = read_param_cfg(orbit_config_struct,param_keyword);
test_georb_mode = strcmp(georb_mode,'orbit_mission');
if test_georb_mode == 1
    param_keyword = 'orbiting_objects_mission';
    [orbiting_objects_mission] = read_param_cfg(orbit_config_struct,param_keyword);
    mission_name = orbiting_objects_mission;
else
    param_keyword = 'orbiting_object_name';
    [orbiting_object_name] = read_param_cfg(orbit_config_struct,param_keyword);
    satellite_id_name = orbiting_object_name;
    test_satellite_mission = strcmp(satellite_id_name,'GRACE-C');
    if test_satellite_mission == 1
       mission_name = 'GRACE_FO_mission'; 
    end
    test_satellite_mission = strcmp(satellite_id_name,'GRACE-D');
    if test_satellite_mission == 1
       mission_name = 'GRACE_FO_mission'; 
    end
    test_satellite_mission = strcmp(satellite_id_name,'GRACE-A');
    if test_satellite_mission == 1
       mission_name = 'GRACE_mission'; 
    end
    test_satellite_mission = strcmp(satellite_id_name,'GRACE-B');
    if test_satellite_mission == 1
       mission_name = 'GRACE_mission'; 
    end    
end

% Satellite missions cases :: GRACE missions
test_mission_grace   = strcmp(mission_name,'GRACE_mission');
test_mission_gracefo = strcmp(mission_name,'GRACE_FO_mission');
% if test_mission_grace == 1 || test_mission_gracefo == 1

% data_matrix = orbc;
% data_functional = 'Orbit';
% reference_frame = 'ICRF';
% [georb2gnv_data_name] = write_georb_data3_gnv(orbit_config_struct, data_matrix, data_functional, reference_frame);
% [status,message,messageid] = movefile(georb2gnv_data_name, OUT_fname_object_mjd);

data_matrix = orbt;
data_functional = 'Orbit';
reference_frame = 'ITRF';
[georb2gnv_data_name] = write_georb_data3_gnv(orbit_config_struct, data_matrix, data_functional, reference_frame);
[status,message,messageid] = movefile(georb2gnv_data_name, OUT_fname_object_mjd);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
