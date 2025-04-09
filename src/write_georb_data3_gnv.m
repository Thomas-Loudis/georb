function [georb2gnv_data_name] = write_georb_data3_gnv(orbit_config_fname, data_matrix, data_functional, reference_frame)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  write_georb_data2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Write computed data to output files in GEORB format e.g. orbits, partial
%  derivatives, observation residuals, external orbit comparison
%  differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - orbit_config_fname  : GEORB configuration structure name 
%
% Output arguments:
% - georb_data_name     : GEORB written data file name
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                               9 July  2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 09/11/2022  Dr. Thomas Loudis Papanikolaou
%             Minor update for introducing write_data_name function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output data file name conventions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[georb_dataformat_name, georb_dataformat_suffix] = write_data_name(orbit_config_fname, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data functional name part1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data_functional_part1 = 'orbit';
data_functional_part1 = data_functional;

% Orbit
data_functional_keyword = 'Orbit';
test_data_functional = strcmp(data_functional,data_functional_keyword);
if test_data_functional == 1
data_functional_part1 = 'orbit';
end

% Partial Derivatives
data_functional_keyword = 'Partial Derivatives';
test_data_functional = strcmp(data_functional,data_functional_keyword);
if test_data_functional == 1
data_functional_part1 = 'partials';
end

% Orbit observation residuals
data_functional_keyword = 'Orbit Observations residuals';
test_data_functional = strcmp(data_functional,data_functional_keyword);
if test_data_functional == 1
data_functional_part1 = 'obs_residuals';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data functional name part2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_functional_part2 = 'crf';

% ICRF
reference_frame_keyword = 'ICRF';
test_reference_frame = strcmp(reference_frame,reference_frame_keyword);
if test_reference_frame == 1
data_functional_part2 = 'crf';
end

% ITRF
reference_frame_keyword = 'ITRF';
test_reference_frame = strcmp(reference_frame,reference_frame_keyword);
if test_reference_frame == 1
data_functional_part2 = 'trf';
end

% ITRF
reference_frame_keyword = 'Kepler';
test_reference_frame = strcmp(reference_frame,reference_frame_keyword);
if test_reference_frame == 1
data_functional_part2 = 'Kepler';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data functional file name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_functional_filename = sprintf('%s%s%s%s','_',data_functional_part1,'_',data_functional_part2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case orbit in itrf
data_functional_keyword = '_orbit_trf';
test_data_functional = strcmp(data_functional_filename,data_functional_keyword);
if test_data_functional == 1
    data_functional_filename = '';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External Orbit Comparison differences 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_functional_keyword = 'External Orbit comparison';
test_data_functional = strcmp(data_functional,data_functional_keyword);
if test_data_functional == 1
    % Orbital Frame
    reference_frame_keyword = 'Orbital Frame';
    test_reference_frame = strcmp(reference_frame,reference_frame_keyword);
    if test_reference_frame == 1
        data_functional_filename = '_delta-RTN';
    end
    % ICRF
    reference_frame_keyword = 'ICRF';
    test_reference_frame = strcmp(reference_frame,reference_frame_keyword);
    if test_reference_frame == 1
        data_functional_filename = '_delta-CRF';
    end
    % Kepler
    reference_frame_keyword = 'Kepler';
    test_reference_frame = strcmp(reference_frame,reference_frame_keyword);
    if test_reference_frame == 1
        data_functional_filename = '_delta-Kepler';
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-satellite observation residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_functional_keyword = 'LRI range-rate residuals';
test_data_functional = strcmp(data_functional,data_functional_keyword);
if test_data_functional == 1
data_functional_filename = '_LRI_rangerate_residuals';
end

data_functional_keyword = 'LRI range residuals';
test_data_functional = strcmp(data_functional,data_functional_keyword);
if test_data_functional == 1
data_functional_filename = '_LRI_range_residuals';
end

data_functional_keyword = 'KBR range-rate residuals';
test_data_functional = strcmp(data_functional,data_functional_keyword);
if test_data_functional == 1
data_functional_filename = '_KBR_rangerate_residuals';
end

data_functional_keyword = 'KBR range residuals';
test_data_functional = strcmp(data_functional,data_functional_keyword);
if test_data_functional == 1
data_functional_filename = '_KBR_range_residuals';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data file name for GNV format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% georb_data_name = sprintf('%s%s%s',georb_dataformat_name,data_functional_filename,georb_dataformat_suffix);

% GNP1B_2019-01-01_C_00.txt

% GEORB to GNV data format name
georb2gnv_dataformat_name = sprintf('%s','GNP1B');

% Date
% mjd_to = data_matrix(1,1);
% Read cofig file for MJDo
[orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(orbit_config_fname);
% fix(IC_MJDo)
mjd_to = IC_MJDo;
[UT,day,month,year] = MJD_inv(mjd_to);
% param_keyword = 'Date';
% fprintf(fid,'%-33s %s %d %d %d ',param_keyword, ': ',day,month,year);
% fprintf(fid,'%s\n','');
% [orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(cfg_filename);
if month < 10
    month_2d = sprintf('%s%d','0',month);
else
    month_2d = sprintf('%d',month);
end
if day < 10
    day_2d = sprintf('%s%d','0',day);
else
    day_2d = sprintf('%d',day);
end
date_name = sprintf('%d%s%s%s%s',year,'-',month_2d,'-',day_2d);

% GRACE/GRACE-FO satellite id
% Satellite/Object name
param_keyword = 'orbiting_object_name';
[orbiting_object_name] = read_param_cfg(orbit_config_fname,param_keyword);
test = strcmp(orbiting_object_name,'GRACE-A');
if test == 1
    grace_sat_id = "A";
end
test = strcmp(orbiting_object_name,'GRACE-B');
if test == 1
    grace_sat_id = "B";
end
test = strcmp(orbiting_object_name,'GRACE-C');
if test == 1
    grace_sat_id = "C";
end
test = strcmp(orbiting_object_name,'GRACE-D');
if test == 1
    grace_sat_id = "D";
end

% Data format suffix
dataformat_suffix = sprintf('%s','.txt');
% dataformat_suffix = georb_dataformat_suffix;

% Data file name in GNV format
georb2gnv_data_name = sprintf('%s', georb2gnv_dataformat_name,'_',date_name,'_',grace_sat_id,'_','01',data_functional_filename,dataformat_suffix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB version
param_keyword = 'GEORB_version';
[src_version] = read_param_cfg(orbit_config_fname, param_keyword);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Write :: data matrix in georb format
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [fid] = write_georb_format(georb_data_name, orbit_config_fname, data_matrix, data_functional, src_version, reference_frame);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write :: data matrix in gnv format (GRACE missions data format)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fid] = write_georb_format_gnv1b(georb2gnv_data_name, orbit_config_fname, data_matrix, data_functional, src_version, reference_frame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
