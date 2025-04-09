function [georb_data_name] = write_georb_statistics(orbit_config_fname, orbit_config_fname_pair, rms_obs, rms_extorb, rms_kbr, rms_lri)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  write_georb_statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Write orbit statistics to data files 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - out_filename  : Orbit output file name
% - orbit_config_fname  : Configuration file name / structure
%
% Output arguments:
% -     : 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                               9 July  2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5/11/2022  Dr. Thomas Loudis Papanikolaou
%            Minor changes due to upgrade of the orbit configuration format based on structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mode: georb_mode :: 'orbit_mission' : Orbit determination of a single orbiting object (satellite, invidual orbiter)
param_keyword = 'georb_mode';
[georb_mode] = read_param_cfg(orbit_config_fname,param_keyword);
test_orbit_mission = strcmp(georb_mode,'orbit_mission');

% Name ID of "satellite mission", "satellite" or "orbiting object"  
param_keyword = 'orbiting_objects_mission';
[orbiting_object_mission_name] = read_param_cfg(orbit_config_fname,param_keyword);

% Satellite missions cases ::
test_grace   = strcmp(orbiting_object_mission_name,'GRACE_mission');
test_gracefo = strcmp(orbiting_object_mission_name,'GRACE_FO_mission');

% Satellite/Object name
param_keyword = 'orbiting_object_name';
[orbiting_object_name] = read_param_cfg(orbit_config_fname,param_keyword);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read cofig file for MJDo
[orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(orbit_config_fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output file name for writing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mission or Object statistics
test_satpair = strcmp(orbit_config_fname_pair,'no_pair');
mission_01 = 0;
if test_satpair == 0
    mission_01 = 1;
end    
[georb_dataformat_name, georb_dataformat_suffix] = write_data_name(orbit_config_fname, mission_01);
georb_dataformat_suffix = sprintf('%s','.out');
georb_data_name = sprintf('%s%s%s',georb_dataformat_name,'_statistics',georb_dataformat_suffix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open file for writing
out_filename = georb_data_name;
fid = fopen(out_filename,'w');

param_keyword = 'GEORB Orbit Statistics';
fprintf(fid,'%-31s',param_keyword);
fprintf(fid,'%s\n\n','');

fprintf(fid,'%-31s \n','-------------------------');

fprintf(fid,'%-31s','RMS');
fprintf(fid,'%s\n','');

fprintf(fid,'%-31s \n','-------------------------');
%fprintf(fid,'%s\n','');

% Format of residuals' text line
format_rms_orb = '%-30s%s %-21s %s %-11.4f %-11.4f %-11.4f';
format_rms_lri = '%-30s%s %-21s %s %-10.6f';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite/Object name read
fprintf(fid,'%-30s%s %-17s','Satellite/Object',':', orbiting_object_name);
fprintf(fid,'%s\n','');

fprintf(fid,format_rms_orb,'Pseudo-Observations residuals',':','ICRF XYZ (m)',':',rms_obs(1,1), rms_obs(1,2), rms_obs(1,3) );
fprintf(fid,'%s\n','');

fprintf(fid,format_rms_orb,'External orbit comparison',':','Orbital Frame RTN (m)',':',rms_extorb(1,1), rms_extorb(1,2), rms_extorb(1,3) );
fprintf(fid,'%s\n\n','');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if test_satpair == 0
% GRACE paired satellite statistcs
if test_grace == 1 || test_gracefo == 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite/Object name read
param_keyword = 'orbiting_object_name';
[orbiting_object_name] = read_param_cfg(orbit_config_fname_pair,param_keyword);

fprintf(fid,'%-30s%s %-17s','Satellite/Object',':', orbiting_object_name);
fprintf(fid,'%s\n','');

fprintf(fid,format_rms_orb,'Pseudo-Observations residuals',':','ICRF XYZ (m)',':',rms_obs(2,1), rms_obs(2,2), rms_obs(2,3) );
fprintf(fid,'%s\n','');

if length(rms_extorb) > 1
fprintf(fid,format_rms_orb,'External orbit comparison',':','Orbital Frame XYZ (m)',':',rms_extorb(2,1), rms_extorb(2,2), rms_extorb(2,3) );
fprintf(fid,'%s\n\n','');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if test_grace == 1 
    fprintf(fid,'%-31s \n','GRACE inter-satellite data residuals');
elseif test_gracefo == 1
    fprintf(fid,'%-31s \n','GRACE-FO inter-satellite data residuals');    
end

% GRACE Intersatellite observations statistcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-band ranging residuals RMS
if length(rms_kbr) > 1
    % Convert to μm/sec
    kbr_rangerate = rms_kbr(1,2) * 10^6; 
    fprintf(fid,format_rms_lri,'KBR data residuals',':','range-rate (μm/sec)',':', kbr_rangerate);
    fprintf(fid,'%s\n','');
end

% LRI ranging residuals RMS
if length(rms_lri) > 1
    % Convert to μm/sec
    lri_rangerate = rms_lri(1,2) * 10^6;  
    fprintf(fid,format_rms_lri,'LRI data residuals',':','range-rate (μm/sec)',':', lri_rangerate);
    fprintf(fid,'%s\n','');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
