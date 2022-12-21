function [orbarray] = sp3c_orb(filename,orbtype,tstop)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function : sp3c_orb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Orbits SP3c format  (e.g. GOCE PSO Kin or Reduced-Dynamic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - filename : Orbit data file name
% - orbtype  : Kinematic or Reduced-Dynamic orbit
%              'kin' or 'dyn'
%
% Output arguments:
% - orbarray:  z = [t r']  or  z = [t r' v']  at every epoch
%   position:  r = [x y z]'       in m
%   velocity:  v = [Vx Vy Vz]'    in m/sec
%   epoch:     t = MJD in Terrestrial Time (TT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         July 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Reading or not the whole day data
if tstop == 0
    tstop = 24 * 3600;
end


fid = fopen(filename);
orbarray_ith = 0;
lp_tstop = 1;
while (~feof(fid)) && lp_tstop == 1
    ln_ith = fgetl(fid);
    if ln_ith(1,1) == '*'
        % GPS time scale
        Year  = sscanf(ln_ith,'%*s %d %*');
        month = sscanf(ln_ith,'%*s %*d %d %*');
        day   = sscanf(ln_ith,'%*s %*d %*d %d %*');
        hour  = sscanf(ln_ith,'%*s %*d %*d %*d %d %*');
        min   = sscanf(ln_ith,'%*s %*d %*d %*d %*d %d %*');
        sec_i = sscanf(ln_ith,'%*s %*d %*d %*d %*d %*d %f %*');
        % sec since 0h
        sec = sec_i + min * 60 + hour * 3600;
        
        % TT time scale & MJD
        [jd_gps,mjd_gps] = MJD_date(TT,day,month,Year);
        [UTC,TT] = time_scales_GPS(sec,mjd_gps);
        [jdTT,mjdTT] = MJD_date(TT,day,month,Year);
        clear jdTT
        
        if TT > tstop 
            lp_tstop = 0;
        end
            
        if UTC >  0
            % Position vector
            ln_ith = fgetl(fid);
            r_ith = sscanf(ln_ith,'%*s %f %f %f');
            r_ith = 10^3 * r_ith;

            if orbtype == 'dyn'
                % Velocity vector
                ln_ith = fgetl(fid);
                v_ith = sscanf(ln_ith,'%*s %f %f %f');
                v_ith = v_ith / 10;
                % Orbit output array
                orbarray_ith = orbarray_ith + 1;
                orbarray(orbarray_ith,:) = [mjdTT r_ith' v_ith'];
            elseif orbtype == 'kin'
                % Orbit output array
                orbarray_ith = orbarray_ith + 1;
                orbarray(orbarray_ith,:) = [mjdTT r_ith'];            
            end       
        end
        
    end
end
fclose(fid);
clear fid
