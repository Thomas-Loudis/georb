function [OUT_fname_object_mjd, OUT_fname_mission_mjd] = write_results_dir(orbit_config_fname,orbit_model_struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: write_results_dir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  write_results_dir forms the output directory of the results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - config_filename     : Orbit Configuration file name
% - results_foldername  : Folder name of the results in-progress being
%                         written prior moving to the final directory path
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             23 August 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Object Directory: naming based on satellite/object name and MJD of initial epoch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name ID of "satellite mission", "satellite" or "orbiting object"  
param_keyword = 'orbiting_object_name';
[orbiting_object_name] = read_param_cfg(orbit_config_fname,param_keyword);

% Read config file for MJDo
IC_MJDo = orbit_model_struct.IC_MJD;

% Directory name according to the satellite/object name and MJDo
OUT_fname_object_mjd = sprintf('%s%s%d', orbiting_object_name,'_',fix(IC_MJDo));

% % Make directory
% [status, message, messageid] = rmdir(OUT_fname_object_mjd,'s');
% [status, message, messageid] = mkdir(OUT_fname_object_mjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move written GEORB output files to folder OUT_fname
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move files to directory
% [status,message,messageid] = movefile('*.orb',OUT_fname);
% [status,message,messageid] = movefile('*.out',OUT_fname);
% [status,message,messageid] = movefile('*.emp',OUT_fname);

% Copy configuration files to directory
%[status,message,messageid] = copyfile(orbit_config_fname ,OUT_fname);
%delete(orbit_config_fname);

% [status,message,messageid] = movefile(data_name, OUT_fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mission Output Directory 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB mode 
param_keyword = 'georb_mode';
[georb_mode] = read_param_cfg(orbit_config_fname,param_keyword);
test_orbit_mission = strcmp(georb_mode,'orbit_mission');

% Objects/Mission :: orbiting_objects_mission
param_keyword = 'orbiting_objects_mission';
[orbiting_objects_mission] = read_param_cfg(orbit_config_fname,param_keyword);

% Satellite missions cases ::
test_mission_grace   = strcmp(orbiting_objects_mission,'GRACE_mission');
test_mission_gracefo = strcmp(orbiting_objects_mission,'GRACE_FO_mission');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE folder for output files/folders of orbits and instersatellite-ranging data analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read config file for MJDo
%[orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(orbit_config_G1);
if test_mission_grace == 1
    %OUT_foldername_GRACE = 'GRACE';
    OUT_fname_mission_mjd = sprintf('%s%s%d','GRACE','_',fix(IC_MJDo));
elseif test_mission_gracefo == 1
    %OUT_foldername_GRACE = 'GRACE-FO';
    OUT_fname_mission_mjd = sprintf('%s%s%d','GRACE-FO','_',fix(IC_MJDo));
end

if test_orbit_mission == 0
    OUT_fname_mission_mjd = 'results_in-progress';     
end    
if test_orbit_mission == 1 
    % [status, message, messageid] = rmdir(OUT_fname_mission_mjd);
    % [status, message, messageid] = mkdir(OUT_fname_mission_mjd);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OUT_fname_results_inprogress = results_foldername;
% %OUT_fname_results_inprogress = 'results_in-progress';
% 
% [status, message, messageid] = rmdir(OUT_fname_results_inprogress);
% [status, message, messageid] = mkdir(OUT_fname_results_inprogress);
% [status,message,messageid] = movefile('*.out',OUT_fname_results_inprogress,'f');
% [status,message,messageid] = movefile(OUT_fname,OUT_fname_results_inprogress,'f');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


