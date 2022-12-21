function [rms_orbital,rms_orbc,rms_orbk,rms_orbt,dstn,dorbc,dkepl,dorbt,rms_3D] = mainf_statistout2(GM,orbc,orbk,orbt,orbce,orbte,veqZarray,veqParray,MODEid,itr_id)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Statistics / Residuals / Orbit Comparison
%  Write statistical quantities to output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments
% - orbc:   State vector in the Celestial Reference System GCRS
%   orbc = [ t r_GCRS' v_GCRS' er' ev']
%   t:      Epoch in seconds of TT (Terrestrial Time)
% - 
% - MODEid : type of orbit analysis
%   >>       EQM : Propagator
%            VEQ : VEQ integration
%            ESM : Estimator
% - itr_id : number of iterations of the estimator procedure
% - fixmjd : Epoch conversion (1st column)  "MJD" to "MJD and t (sec)"
%            MJD is including the fraction of the day
%            MJD convert or not to integer
%   fixmjd = 1 : convert MJD to integer
%   fixmjd = 0 : do not convert MJD to integer

%
% Output arguments:
% - orbc:         State vector in the Celestial Reference System GCRS
%   orbc = [ t r_GCRS' v_GCRS' er' ev']
%   t:            Epoch in seconds of TT (Terrestrial Time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                           May 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Width number
wno = 29;
% Precision : number of decimal digits
pr_orbc = [0 9 11 11 11 14 14 14]';
pr_orbk = [0 9 9 12 9 9 9 9]';
pr_RESorbital = [0 9 9 9 9]';
pr_rmsorbc = [9 9 9 12 12 12]';
pr_rmsorbk = [9 10 8 8 8 8]';
pr_rmsorbital = [9 9 9]';
pr_veqZ = [0 9 9 9 9 9 9 9]';
pr_veqP = [9 9 9 9 9 9 9 9]';
% Epoch conversion (1st column) : "MJD" to "MJD and t (sec)"
fixmjd = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for Velocity vector
[d1 d2] = size(orbte);
if d2 < 7
    Vel_vector = 0;
elseif d2 >= 7
    Vel_vector = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Vel_vector == 1
    % Orbital perturbations and differences
    [dstn,rms_orbital,dorbc,rms_orbc,dkepl,rms_orbk] = orbital_pert(orbce,orbc,GM);
    % Terrestrial Reference frame (ITRS)
    [dorbt,rms_orbt] = compstat(orbte,orbt);
    % set 3D RMS to 0 for Vel_vector=1
    rms_3D = 0;
elseif Vel_vector == 0
    orbt_pos = orbt(:,1:4);
    orbc_pos = orbc(:,1:4);
    % Terrestrial Reference frame (ITRS)
    [dorbt,rms_orbt,orbt_common] = compstat(orbte,orbt_pos);
    % Celestial Reference Frame (GCRS)
    [dorbc,rms_orbc,orbc_common] = compstat(orbce,orbc_pos);  
    % 3D RMS (ITRF and ICRF)
    % proved equation    rms_3d = sqrt(rms_x^2 + rms_y^2 + rms_z^2)
    RMS_ITRF_3D = sqrt(rms_orbt(1,1)^2 + rms_orbt(1,2)^2 + rms_orbt(1,3)^2);
    RMS_ICRF_3D = sqrt(rms_orbc(1,1)^2 + rms_orbc(1,2)^2 + rms_orbc(1,3)^2);
    rms_3D = [RMS_ITRF_3D RMS_ICRF_3D];
    
    % Set Zeros for other Frames
    dkepl = 0;
    rms_orbk = 0;
    dstn = 0;
    rms_orbital = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write data to files condition:
%if itr_id > -1
if -2 > -1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
MODEid2write = sprintf('%s%d',MODEid,itr_id);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epoch conversion (1st column) : "MJD" to "MJD and t (sec)"
[orbc2] = mjd2mjdtt(orbc,fixmjd);
[orbt2] = mjd2mjdtt(orbt,fixmjd);
[orbk2] = mjd2mjdtt(orbk,fixmjd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write data to files: Orbit Determination output data 
outfilename1 = sprintf('%s%s%s','orbk',MODEid2write,'.out');
outfilename2 = sprintf('%s%s%s','orbc',MODEid2write,'.out');
outfilename3 = sprintf('%s%s%s','orbt',MODEid2write,'.out');
[formatarray] = writedat(outfilename1,wno,pr_orbk,orbk2');
[formatarray] = writedat(outfilename2,wno,pr_orbc,orbc2');
[formatarray] = writedat(outfilename3,wno,pr_orbc,orbt2');
clear outfilename1 outfilename2 outfilename3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write data to files: external orbit data
[orbce_2] = mjd2mjdtt(orbce,fixmjd);
[orbte_2] = mjd2mjdtt(orbte,fixmjd);
outfilename2 = sprintf('%s%s','orbc_ext','.out');
outfilename3 = sprintf('%s%s','orbt_ext','.out');
[formatarray] = writedat(outfilename2,wno,pr_orbc,orbce_2');
[formatarray] = writedat(outfilename3,wno,pr_orbc,orbte_2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write data to files: RMS report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outfilename = sprintf('%s%s%s','STATreport_',MODEid2write,'.out');
fid = fopen(outfilename,'w');
fprintf(fid,'%s\n','-------------------------');
fprintf(fid,'%s\n','Orbit Statistics');
fprintf(fid,'%s\n\n','-------------------------');
fprintf(fid,'%s\n','RMS');
fprintf(fid,'%s\n','-------------------------');
if Vel_vector == 1
    [formatarray0] = writeformat(wno,pr_rmsorbc);
    formatarray = sprintf('%s %s','%s',formatarray0);
    fprintf(fid,formatarray,'ICRF                    :',rms_orbc');
    fprintf(fid,formatarray,'ITRF                    :',rms_orbt');
    [formatarray0] = writeformat(wno,pr_rmsorbk);
    formatarray = sprintf('%s %s','%s',formatarray0);
    fprintf(fid,formatarray,'Kepler Elements         :',rms_orbk');
    [formatarray0] = writeformat(wno,pr_rmsorbital);
    formatarray = sprintf('%s %s','%s',formatarray0);
    fprintf(fid,formatarray,'Orbital Frame           :',rms_orbital');
elseif Vel_vector == 0
    pr_rmsXYZ = pr_rmsorbc(1:3,1);
    [formatarray0] = writeformat(wno,pr_rmsXYZ);
    formatarray = sprintf('%s %s','%s',formatarray0);
    fprintf(fid,formatarray,'ICRF                    :',rms_orbc');
    fprintf(fid,formatarray,'ITRF                    :',rms_orbt');
%     fprintf(fid,'%s %17.6f %17.6f','3D                      :',rms_3D);
    fprintf(fid,'%s %17.9f \n','ICRF 3D                 :',RMS_ICRF_3D);
    fprintf(fid,'%s %17.9f \n','ITRF 3D                 :',RMS_ITRF_3D);
end
fclose(fid);
clear outfilename formatarray0 formatarray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Residuals / Orbit comparison
% Epoch conversion (1st column) : "MJD" to "MJD and t (sec)"
outfilename1 = sprintf('%s%s%s','RESorbk_',MODEid2write,'.out');
outfilename2 = sprintf('%s%s%s','RESorbc_',MODEid2write,'.out');
outfilename3 = sprintf('%s%s%s','RESorbt_',MODEid2write,'.out');
outfilename4 = sprintf('%s%s%s','RESorbital_',MODEid2write,'.out');

if Vel_vector == 0
    pr_orbc = [0 9 4 4 4]';
end
[dorbc2] = mjd2mjdtt(dorbc,fixmjd);
[dorbt2] = mjd2mjdtt(dorbt,fixmjd);
[formatarray] = writedat(outfilename2,wno,pr_orbc,dorbc2');
[formatarray] = writedat(outfilename3,wno,pr_orbc,dorbt2');

if Vel_vector == 1
    [dstn2] = mjd2mjdtt(dstn,fixmjd);
    [dorbk2] = mjd2mjdtt(dkepl,fixmjd);
    [formatarray] = writedat(outfilename1,wno,pr_orbk,dorbk2');
    [formatarray] = writedat(outfilename4,wno,pr_RESorbital,dstn2');
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write data to files: VEQ arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz sz2] = size(veqZarray);
if sz > 1
    [veqZarray2] = mjd2mjdtt(veqZarray,fixmjd);
    outfilename1 = sprintf('%s%s%s','veqZarray',MODEid2write,'.out');
    [formatarray] = writedat(outfilename1,wno,pr_veqZ,veqZarray2');
    clear outfilename1 outfilename2
end

% VEQ P matrix
    [veqParray2] = mjd2mjdtt(veqParray,fixmjd);
    outfilename1 = sprintf('%s%s%s','veqParray',MODEid2write,'.out');
    [formatarray] = writedat(outfilename1,wno,pr_veqP,veqParray2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
