function georb_series(orbit_model_struct, ic_data_struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: georb_series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  georb_series performs series of orbit determination for all objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - orbit_model_struct  : Orbit model structure array 
% - ic_data_struct      : Initial Conditions structure array 
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                            21 January 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Series of georb_function computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ic_struct_d1, ic_struct_d2] = size(ic_data_struct);
[ic_n, ic_m] = size(ic_data_struct(1).ic_data);

for ic_i = 1 : ic_n
    for constellation_id = 1 : ic_struct_d2
        ic_data_struct_run(constellation_id).ic_data = ic_data_struct(constellation_id).ic_data(ic_i,:);
    end
    [out_dir_name] = georb_function(orbit_model_struct, ic_data_struct_run);
end
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
