function [georb_data_name] = write_georb_data2(orbit_config_fname, data_matrix, data_functional, reference_frame)


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
data_functional_part1 = 'orbit';

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
% Data file name
georb_data_name = sprintf('%s%s%s',georb_dataformat_name,data_functional_filename,georb_dataformat_suffix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB version
param_keyword = 'GEORB_version';
[src_version] = read_param_cfg(orbit_config_fname, param_keyword);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write :: data matrix in georb format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fid] = write_georb_format(georb_data_name, orbit_config_fname, data_matrix, data_functional, src_version, reference_frame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

