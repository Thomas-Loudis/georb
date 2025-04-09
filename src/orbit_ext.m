function [rms_orbital,rms_orbc,rms_orbk,rms_orbt,dstn,dorbc,dkepl,dorbt,orbce,orbte] = orbit_ext(orbc,orbk,orbt,veqZarray,veqParray,MODEid,cfg_fname, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_ext.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  External orbit comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname :  Orbit oncfiguration file (format 2022) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         June 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 08/07/2022   Thomas Loudis Papanikolaou 
%              Code upgrade based on the new orbit configuration format
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit/Force model central matrix
GM_glob = orbit_model_struct.GM_Earth;
eopdat  = orbit_model_struct.EOP_data;
dpint   = orbit_model_struct.EOP_interp_no;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write to output files y/n :: 1/0
write_out_flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf parameters:
% Width number
wno = 29;
% Precision : number of decimal digits
pr_orbc = [0 9 11 11 11 14 14 14]';
pr_orbk = [0 9 9 12 9 9 9 9]';
pr_rmsorbc = [0 9 9 9 12 12 12]';
pr_rmsorbk = [0 9 12 9 9 9 9]';
pr_rmsorbital = [0 9 9 9]';
pr_veqZ = [0 9 9 9 9 9 9 9]';
% Epoch conversion (1st column) : "MJD" to "MJD and t (sec)"
fixmjd = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read configuration file for the option of external orbit comparison
param_keyword = 'external_orbit_comp';
[param_value] = read_param_cfg(cfg_fname,param_keyword);
extorb_test = param_value;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External Orbit 
if extorb_test == 'y'
    % Read external orbit data
    [orbte,orbce,orbke, ext_orbit_comp_yn] = prm_orbext(cfg_fname,GM_glob,eopdat,dpint,orbit_model_struct);
    if write_out_flag > 0
    [orbke2] = mjd2mjdtt(orbke,fixmjd);
    [orbce2] = mjd2mjdtt(orbce,fixmjd);
    [orbte2] = mjd2mjdtt(orbte,fixmjd);
    % Write data to files: orbits (external)
    outfilename1 = 'orbke.out';
    outfilename2 = 'orbce.out';
    outfilename3 = 'orbte.out';
    [formatarray] = writedat(outfilename1,wno,pr_orbk,orbke2');
    [formatarray] = writedat(outfilename2,wno,pr_orbc,orbce2');
    [formatarray] = writedat(outfilename3,wno,pr_orbc,orbte2');
    end
elseif extorb_test == 'n'
    rms_orbital = zeros(1,3);    
    rms_orbc = zeros(1,6);
    rms_orbk = rms_orbc;
    rms_orbt = rms_orbc;
    dstn  = zeros(1,4);
    dorbc = zeros(1,7);
    dkepl = zeros(1,7);
    dorbt = zeros(1,7);
    orbce = zeros(1,7);
    orbte = zeros(1,7);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit comparison between computed orbit and external orbit
if extorb_test == 'y'
    % External Orbit Comparison / Write Data
    MODEid2 = sprintf('%s%s',MODEid,'extorb');
    i_itr = -1;
    [rms_orbital,rms_orbc,rms_orbk,rms_orbt,dstn,dorbc,dkepl,dorbt,rms_3D] = mainf_statistout2(GM_glob,orbc,orbk,orbt,orbce,orbte,veqZarray,veqParray,MODEid2,i_itr);
elseif extorb_test == 'n'
    if write_out_flag > 0
        % Write data to files: orbits [orbc orbk orbt]
        %MODEid2write = sprintf('%s%d',MODEid,i_itr);
        MODEid2write = sprintf('%s',MODEid);
        % Epoch conversion (1st column) : "MJD" to "MJD and t (sec)"
        [orbk2] = mjd2mjdtt(orbk,fixmjd);
        [orbc2] = mjd2mjdtt(orbc,fixmjd);
        [orbt2] = mjd2mjdtt(orbt,fixmjd);
        outfilename1 = sprintf('%s%s%s','orbk',MODEid2write,'.out');
        outfilename2 = sprintf('%s%s%s','orbc',MODEid2write,'.out');
        outfilename3 = sprintf('%s%s%s','orbt',MODEid2write,'.out');
        [formatarray] = writedat(outfilename1,wno,pr_orbk,orbk2');
        [formatarray] = writedat(outfilename2,wno,pr_orbc,orbc2');
        [formatarray] = writedat(outfilename3,wno,pr_orbc,orbt2');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


