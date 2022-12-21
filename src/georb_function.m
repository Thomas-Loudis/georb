function georb_function(main_config_fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: georb_mainfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  georb_mainfunction is the main function for calling the source code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - config_filename     : Configuration file name
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                               5 July  2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 05/07/2022  Thomas Loudis Papanikolaou
%             Code modification for changing the script to function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


src_version = 'v.1.7.5';
fprintf('%s%s \n\n','GEORB ',src_version);

to_tic = tic;
delete('*.out');
[status, message, messageid] = rmdir('OUT*','s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB main settings :: Read main config file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB mode 
param_keyword = 'georb_mode';
[georb_mode] = read_param_file(main_config_fname,param_keyword);

% Name ID of "satellite mission", "satellite" or "orbiting object"  
param_keyword = 'orbiting_object_name';
[orbiting_object_name] = read_param_file(main_config_fname,param_keyword);

% Orbit modelling configuration files
param_keyword = 'orb_config_filename';
[orbit_model_filename] = read_param_file(main_config_fname,param_keyword);

% IC configuration file
param_keyword = 'ic_config_filename';
[ic_config_filename] = read_param_file(main_config_fname,param_keyword);

% % Satellite Data configuration file
% param_keyword = 'sat_data_filename';
% [sat_data_filename] = read_param_file(main_config_fname,param_keyword);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove previously written files and folders
delete('*.out');
delete('*.orb');
delete('*.emp');
[status, message, messageid] = rmdir('results*','s');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precise Orbit Determination 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Mode :: Satellites/Objects series orbit determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Mode: georb_mode :: 'orbit_objects' : Orbit determination of individual orbiting objects (satellite, orbiter)
test_orbit_objects = strcmp(georb_mode,'orbit_objects');
   
if test_orbit_objects == 1   
    % Orbit Determination of objects/satellites series
    [out_dir_name] = orbit_object_series(main_config_fname, src_version);
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Mode: GRACE missions orbit determination and intersatellite data analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Mode: georb_mode :: 'orbit_mission' 
test_orbit_mission = strcmp(georb_mode,'orbit_mission');

% Satellite missions cases ::
test_mission_grace   = strcmp(orbiting_object_name,'GRACE_mission');
test_mission_gracefo = strcmp(orbiting_object_name,'GRACE_FO_mission');

if test_orbit_mission == 1 && ( test_mission_grace == 1 || test_mission_gracefo == 1 )
    [out_dir_name] = orbit_mission_grace(main_config_fname, src_version);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Mode :: Time series of orbit determination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Mode: 'orbit_mission_series' :: Time series of GRACE missions orbit determination
test_orbit_mission_series = strcmp(georb_mode,'orbit_mission_series');
test_time_series = strcmp(georb_mode,'orbit_time_series');
%if test_orbit_mission_series == 1 && ( test_mission_grace == 1 || test_mission_gracefo == 1 )
if test_time_series == 1    
    fprintf('%s \n', 'Orbit Time series mode is not available by the current release');
    %[out_dir_name] = orbit_time_series(main_config_fname, src_version);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode :: Orbit Design, Orbit Simulation of GRACE-like missions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = strcmp(georb_mode,'orbit_design');
if test == 1
    fprintf('%s \n', 'Orbit Design mode is not available by the current release');
    %[orb_crf_1, orb_trf_1, orb_kep_1, orb_crf_2, orb_trf_2, orb_kep_2, Xmatrix, range_matrix_trf, dRTN_matrix_trf,rms_range_trf] = orbit_design_formation (config_filename);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s % .3f \n', 'Overall Computation Time (min)             :', toc(to_tic)/60);
[status,message,messageid] = copyfile(main_config_fname,out_dir_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Results directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output epoch (computer time)
[time_clock] = time_filename();
% Output: Results Directory name
results_dir_name_0 = 'results';
time_name = time_clock;
results_dir_name = sprintf('%s%s%s',results_dir_name_0,'_',time_name);
%[status, message, messageid] = mkdir('../data_output', out_dir_name);
%[status,message,messageid] = mkdir(results_dir_name);
[status,message,messageid] = movefile(out_dir_name,results_dir_name);
[status,message,messageid] = movefile(results_dir_name,'../data_output');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
