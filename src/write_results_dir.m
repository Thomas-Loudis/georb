function write_results_dir(orbit_config_fname, results_foldername)

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
% Directory name according to the satellite/object name and MJD of initial epoch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name ID of "satellite mission", "satellite" or "orbiting object"  
param_keyword = 'orbiting_object_name';
[orbiting_object_name] = read_param_cfg(orbit_config_fname,param_keyword);

% Read config file for MJDo
[orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(orbit_config_fname);

% Directory name according to the satellite/object name and MJDo
OUT_fname = sprintf('%s%s%d', orbiting_object_name,'_',fix(IC_MJDo));

% Make directory
[status, message, messageid] = rmdir(OUT_fname,'s');
[status, message, messageid] = mkdir(OUT_fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move written GEORB output files to folder OUT_fname
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Move files to directory
[status,message,messageid] = movefile('*.orb',OUT_fname);
[status,message,messageid] = movefile('*.out',OUT_fname);
[status,message,messageid] = movefile('*.emp',OUT_fname);

% Copy configuration files to directory
%[status,message,messageid] = copyfile(orbit_config_fname ,OUT_fname);
%delete(orbit_config_fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Folder 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ic_i == 1
%     OUT_fname_results_inprogress = 'results_in-progress';
%     [status, message, messageid] = rmdir(OUT_fname_results_inprogress);
%     [status, message, messageid] = mkdir(OUT_fname_results_inprogress);
% end
% [status,message,messageid] = movefile('*.out',OUT_fname_results_inprogress,'f');
% [status,message,messageid] = movefile(out_dir_name,OUT_fname_results_inprogress,'f');

OUT_fname_results_inprogress = results_foldername;
%OUT_fname_results_inprogress = 'results_in-progress';

[status, message, messageid] = rmdir(OUT_fname_results_inprogress);
[status, message, messageid] = mkdir(OUT_fname_results_inprogress);
[status,message,messageid] = movefile('*.out',OUT_fname_results_inprogress,'f');
[status,message,messageid] = movefile(OUT_fname,OUT_fname_results_inprogress,'f');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


