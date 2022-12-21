function [fid] = write_georb_data(orbit_config_fname, orbc, orbt, orbk, veqZarray, veqParray, obs_residuals, extorbdiff_stn, extorbdiff_orbc, extorbdiff_kepl)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  write_georb_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Write orbits and partial derivatives to data files 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - out_filename  : Orbit output file name
% - orbit_config_fname  : Configuration file name
%
% Output arguments:
% -     : 
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
% GEORB version
param_keyword = 'GEORB_version';
[src_version] = read_param_cfg(orbit_config_fname, param_keyword);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write :: Orbits to files in georb format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbits
data_functional  = 'Orbit';

% Orbit in ICRF
reference_frame  = 'ICRF';
%georb_orbit_name = sprintf('%s%s%d%s%s', orbiting_object_name,'_',fix(IC_MJDo),'_icrf','.orb');
georb_orbit_name = sprintf('%s%s%s',georb_dataformat_name,'_orbit_crf',georb_dataformat_suffix);
[fid] = write_georb_format(georb_orbit_name, orbit_config_fname, orbc, data_functional, src_version, reference_frame);

% Orbit in ITRF
reference_frame  = 'ITRF';
%georb_orbit_name = sprintf('%s%d%s%s', orbiting_object_name,'_',fix(IC_MJDo),'_itrf','.orb');
georb_orbit_name = sprintf('%s%s%s',georb_dataformat_name,'_orbit_trf',georb_dataformat_suffix);
[fid] = write_georb_format(georb_orbit_name, orbit_config_fname, orbt, data_functional, src_version, reference_frame);

% Orbit: Kepler elements
reference_frame  = 'Kepler';
%georb_orbit_name = sprintf('%s%d%s%s', orbiting_object_name,'_',fix(IC_MJDo),'_kepler','.orb');
georb_orbit_name = sprintf('%s%s%s',georb_dataformat_name,'_kepler',georb_dataformat_suffix);
[fid] = write_georb_format(georb_orbit_name, orbit_config_fname, orbk, data_functional, src_version, reference_frame);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write :: Partial derivatives :: Variational Equations' solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n1 n2] = size(veqZarray);
[n3 n4] = size(veqParray);
if n1 > 1 

if n1 == n3
VEQ_matrix = [veqZarray veqParray(:,2:n4)];
end
if n3 <= 1 
VEQ_matrix = [veqZarray];
end    

% Data functional name
data_functional  = 'Partial Derivatives';
% Data functional file name part
data_functional_filename = '_partials';
% Reference frame
reference_frame  = 'ICRF';
% Data file name
georb_orbit_name = sprintf('%s%s%s',georb_dataformat_name,data_functional_filename,georb_dataformat_suffix);
% Write data in georb format
[fid] = write_georb_format(georb_orbit_name, orbit_config_fname, VEQ_matrix, data_functional, src_version, reference_frame);

% reference_frame  = 'ICRF';
% georb_orbit_name = sprintf('%s%s%s',georb_dataformat_name,'_partials_veqZ',georb_dataformat_suffix);
% [fid] = write_georb_format(georb_orbit_name, orbit_config_fname, veqZarray, data_functional, src_version, reference_frame);
% 
% reference_frame  = 'ICRF';
% georb_orbit_name = sprintf('%s%s%s',georb_dataformat_name,'_partials_veqP',georb_dataformat_suffix);
% [fid] = write_georb_format(georb_orbit_name, orbit_config_fname, veqParray, data_functional, src_version, reference_frame);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write :: Orbit Observation residuals to files in georb format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n1 n2] = size(obs_residuals);
if n1 > 1 
% Data functional name
data_functional  = 'Orbit Observations residuals';
% Data functional file name part
data_functional_filename = '_obs_residuals';
% Orbit residuals in ICRF
reference_frame  = 'ICRF';
% Data file name
georb_orbit_name = sprintf('%s%s%s',georb_dataformat_name,data_functional_filename,georb_dataformat_suffix);
% Write data in georb format
[fid] = write_georb_format(georb_orbit_name, orbit_config_fname, obs_residuals, data_functional, src_version, reference_frame);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extorbdiff_stn, extorbdiff_orbc, extorbdiff_kepl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write :: External Orbit Comparison differences 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_functional  = 'External Orbit comparison';

% Orbital Frame
[n1 n2] = size(extorbdiff_stn);
if n1 > 1 
reference_frame  = 'Orbital Frame';
% Data functional file name part
data_functional_filename = '_delta-RTN';
georb_orbit_name = sprintf('%s%s%s',georb_dataformat_name, data_functional_filename, georb_dataformat_suffix);
[fid] = write_georb_format(georb_orbit_name, orbit_config_fname, extorbdiff_stn, data_functional, src_version, reference_frame);
end

% Orbit differences in ICRF
[n1 n2] = size(extorbdiff_orbc);
if n1 > 1 
reference_frame  = 'ICRF';
% Data functional file name part
data_functional_filename = '_delta-CRF';
georb_orbit_name = sprintf('%s%s%s',georb_dataformat_name,data_functional_filename,georb_dataformat_suffix);
[fid] = write_georb_format(georb_orbit_name, orbit_config_fname, extorbdiff_orbc, data_functional, src_version, reference_frame);
end

% Orbit: Kepler elements
[n1 n2] = size(extorbdiff_kepl);
if n1 > 1 
reference_frame  = 'Kepler';
% Data functional file name part
data_functional_filename = '_delta-Kepler';
georb_orbit_name = sprintf('%s%s%s',georb_dataformat_name,data_functional_filename,georb_dataformat_suffix);
[fid] = write_georb_format(georb_orbit_name, orbit_config_fname, extorbdiff_kepl, data_functional, src_version, reference_frame);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
georb_orbit_name = sprintf('%s%s%s',georb_dataformat_name,'_statistics',georb_dataformat_suffix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
