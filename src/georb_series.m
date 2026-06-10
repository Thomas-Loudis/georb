function georb_series(orbit_model_struct, ic_satellites) 

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
% Function: georb_series 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
%  georb_series performs series of orbit determination for all objects 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - orbit_model_struct  : Orbit model structure array  
% - ic_satellites       : Initial Conditions for all satellites for all days in structure array 
%                       Structure Format: ic_satellites.epochs.ic_data 
% Output arguments: 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                            21 January 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Parallel Computing settings 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
parallel_computing_yn = orbit_model_struct.parallel_computing.parallel_computing_yn; 
num_cpu_cores = orbit_model_struct.parallel_computing.num_cpu_cores; 
 
georb_par_mode = 0; 
test = strcmp(parallel_computing_yn,'y'); 
if test == 1 
    georb_par_mode = 1; 
end 
 
% Number of Cores to be used in parallel 
Ncpus_cores = num_cpu_cores; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Series of georb_function computations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[ic_struct_d1, ic_struct_d2] = size(ic_satellites); 
No_satellites = ic_struct_d2; 
 
ic_epochs = ic_satellites.epochs; 
[ic_epochs_d1, ic_epochs_d2] = size(ic_epochs); 
No_ic_epochs = ic_epochs_d2; 
ic_n = No_ic_epochs; 
 
if georb_par_mode == 0  
 
ic_satellites_run = struct('ic_data',{}); 
 
if 1 < 0 
for ic_i = 1 : ic_n 
    for constellation_id = 1 : No_satellites 
        % ic_satellites_run(constellation_id).ic_data = ic_satellites(constellation_id).ic_data(ic_i,:); 
        ic_satellites_run(constellation_id).ic_data = ic_satellites(constellation_id).epochs(ic_i).ic_data; 
    end 
    % Call main georb_function 
    [out_dir_name] = georb_function(orbit_model_struct, ic_satellites_run); 
end 
end 
 
for ic_i = 1 : ic_n 
% Call georb_par for serial/parallel call of the main function georb_function     
ic_index = ic_i; 
[out_dir_name] = georb_par(orbit_model_struct, ic_satellites, ic_index); 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif georb_par_mode == 1 
parpool("Processes", Ncpus_cores); 
 
% parfor ic_i = 1 : ic_n     
parfor (ic_i = 1 : ic_n, Ncpus_cores)     
% Call georb_par for serial/parallel call of the main function georb_function     
ic_index = ic_i; 
[out_dir_name] = georb_par(orbit_model_struct, ic_satellites, ic_index); 
% georb_par(orbit_model_struct, ic_satellites, ic_index); 
end 
out_dir_name = 'results_in-progress'; 
 
end 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% elseif georb_par_mode == 2 
%  
%  
% % id_parcluster = parcluster() 
% % for ic_i = 1 : ic_n 
% %     for constellation_id = 1 : ic_struct_d2 
% %         % ic_satellites_run(constellation_id).ic_data = ic_satellites(constellation_id).ic_data(ic_i,:); 
% %         ic_satellites_run(constellation_id).ic_data = ic_satellites(constellation_id).epochs(ic_i).ic_data; 
% %     end 
% %     % Call main georb_function 
% %     % [out_dir_name] = georb_function(orbit_model_struct, ic_satellites_run); 
% %     job = batch(id_parcluster,@georb_function,0,{orbit_model_struct, ic_satellites_run}); 
% % end 
% % wait(job); 
% % results_dir_name = 'results_in-progress'; 
% % out_dir_name = results_dir_name; 
%  
%  
% end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Results directory 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Output epoch (computer time) 
[MJD_clock, calendar_clock] = time_filename(); 
% Output: Results Directory name 
results_dir_name_0 = 'results'; 
time_name = calendar_clock; 
results_dir_name = sprintf('%s%s%s',results_dir_name_0,'_',time_name); 
 
% Rename out_dir_name to results_dir_name 
[status,message,messageid] = movefile(out_dir_name,results_dir_name); 
 
% Move written output folders/files to results directory 
results_dir_path = fullfile(pwd,results_dir_name); 
[status,message,messageid] = movefile('*.orb',results_dir_name); 
[status,message,messageid] = movefile('*.out',results_dir_name); 
 
% Data Output main directory 
data_output_path = fullfile(pwd,'..','data_output'); 
data_output_isfolder = isfolder(data_output_path); 
if data_output_isfolder == 0 
[status, message, messageid] = mkdir(data_output_path); 
end 
[status,message,messageid] = movefile(results_dir_name,data_output_path); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
