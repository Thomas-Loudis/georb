function [obsorbc,obsorbt,obsorbc_ext,obsorbt_ext,obsorbc_full,obsorbt_full,COVobs,COVPform] = orbit_obs(cfg_fname)


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables called
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global GM_glob 
global EOP_DAT_glob EOP_dpint_glob
eopdat  = EOP_DAT_glob;
dpint   = EOP_dpint_glob;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observations model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinematic Orbit considered as pseudo-Observations
% OBS : Kinematic Orbit positions
COVPform = 0;
[obsorbc, obsorbk, obsorbt, obsorbc_ext, obsorbk_ext, obsorbt_ext, obsorbc_full, obsorbk_full, obsorbt_full, COVobs] = prm_pseudobs(cfg_fname,GM_glob,eopdat,dpint);
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

