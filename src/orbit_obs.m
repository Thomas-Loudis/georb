function [obsorbc,obsorbt,obsorbc_ext,obsorbt_ext,obsorbc_full,obsorbt_full,COVobs,COVPform, observation_struct] = orbit_obs(cfg_fname, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_obs.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Pseudo-Observations formation based on kinematic orbit data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         June 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 08/07/2022, Thomas Loudis Papanikolaou
%             Code upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orbit model struct matrix 
GM_glob = orbit_model_struct.GM_Earth;
eopdat = orbit_model_struct.EOP_data;
dpint = orbit_model_struct.EOP_interp_no;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observations model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinematic Orbit considered as pseudo-Observations
% OBS : Kinematic Orbit positions
COVPform = 0;
[obsorbc, obsorbk, obsorbt, obsorbc_ext, obsorbk_ext, obsorbt_ext, obsorbc_full, obsorbk_full, obsorbt_full, COVobs] = prm_pseudobs(cfg_fname,GM_glob,eopdat,dpint, orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1 < 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write OBS to files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for Velocity vector
[d1 d2] = size(obsorbc_full);
if d2 < 7
    Vel_vector = 0;
elseif d2 >= 7
    Vel_vector = 1;
end
fixmjd = 1;
if Vel_vector == 0
    [d1 d2] = size(obsorbc);
    obsorbcV = [obsorbc zeros(d1,1) zeros(d1,1) zeros(d1,1)];
    [d1 d2] = size(obsorbc_ext);
    obsorbc_extV = [obsorbc_ext zeros(d1,1) zeros(d1,1) zeros(d1,1)];
    [d1 d2] = size(obsorbc_full);
    obsorbc_fullV = [obsorbc_full zeros(d1,1) zeros(d1,1) zeros(d1,1)];
    [obsorbc2] = mjd2mjdtt(obsorbcV,fixmjd);
    [obsorbc_ext2] = mjd2mjdtt(obsorbc_extV,fixmjd);
    [obsorbc_full2] = mjd2mjdtt(obsorbc_fullV,fixmjd);
elseif Vel_vector == 1
    [obsorbc2] = mjd2mjdtt(obsorbc,fixmjd);
    [obsorbc_ext2] = mjd2mjdtt(obsorbc_ext,fixmjd);
    [obsorbc_full2] = mjd2mjdtt(obsorbc_full,fixmjd);   
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf parameters:
% Width number
wno = 29;
% Precision : number of decimal digits
pr_orbc = [0 9 11 11 11 14 14 14]';
outfilename1 = 'OBS.out';
outfilename2 = 'OBSext.out';
outfilename3 = 'OBSall.out';
[formatarray] = writedat(outfilename1,wno,pr_orbc,obsorbc2');
[formatarray] = writedat(outfilename2,wno,pr_orbc,obsorbc_ext2');
[formatarray] = writedat(outfilename3,wno,pr_orbc,obsorbc_full2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit observations in ICRF 
observation_struct.obs_orbit_crf = obsorbc;
% Orbit observations in ITRF 
observation_struct.obs_orbit_trf = obsorbt;

% Observations' covariance matrix 
observation_struct.obs_cov = COVobs;

% Observations' weight matrix form 
observation_struct.obs_weight = COVPform;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
