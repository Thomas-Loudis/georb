function [orbit_config_struct, orbit_matrix, orbit_rms, veqZarray, veqParray, OBS_matrix, Xparam_aposteriori] = orbit_object(orbit_config_struct, write_data)

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
    fprintf('\n%s\n','Mode: Orbit Determinaton');    
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
fprintf('%s\n','Models and Data preprocessing :: In-progress');    
% Call function orbit_pod
[rms_orbital,rms_orbc,rms_orbk,rms_orbt,dstn,dorbc,dkepl,dorbt,Xmatrix,orbc,orbk,orbt,veqZarray,veqParray,rms_obs,OBS_matrix, Xparam_aposteriori, OBS_residuals] = orbit_pod(orbit_config_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = strcmp(orbit_mode,'orbit_determination');
if test == 1
    fprintf('\n');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbits
data_matrix = orbc;
data_functional = 'Orbit';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);

data_matrix = orbt;
data_functional = 'Orbit';
reference_frame = 'ITRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);

data_matrix = orbk;
data_functional = 'Orbit';
reference_frame = 'Kepler';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial Derivatives
[n1 n2] = size(veqZarray);
[n3 n4] = size(veqParray);
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit observation residuals
[n1 n2] = size(OBS_residuals);
if n1 > 1 
data_matrix = OBS_residuals;
data_functional = 'Orbit Observations residuals';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External Orbit comparison
[n1 n2] = size(dstn);
if n1 > 1 
data_matrix = dstn;
data_functional = 'External Orbit comparison';
reference_frame = 'Orbital Frame';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);

data_matrix = dkepl;
data_functional = 'External Orbit comparison';
reference_frame = 'Kepler';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);

data_matrix = dorbc;
data_functional = 'External Orbit comparison';
reference_frame = 'ICRF';
[georb_data_name] = write_georb_data2(orbit_config_struct, data_matrix, data_functional, reference_frame);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Statistics
rms_kbr = 0;
rms_lri = 0;
orbit_config_fname_pair = 'no_pair';
[fid] = write_georb_statistics(orbit_config_struct, orbit_config_fname_pair, rms_obs, rms_orbital, rms_kbr, rms_lri);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
