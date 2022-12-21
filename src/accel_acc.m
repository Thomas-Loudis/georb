function [facclx,faccly,facclz] = accel_acc(mjd,t,accl,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% accel_acc:  Acceleration by accelerometry data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the satellite acceleration based on accelerometer data
%  The accelerometry represents the non-graviational effects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - mjd:  computation's epoch MJD including fraction of the day
% - t:    computation's time in TT (seconds since 0h)
% - accl: calibrated accelerometry data in inertial system ICRS
%
% Output arguments:
% - Acceleration's cartesian components in ICRS at epoch t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         April 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[sz1 sz2] = size(accl);
%rate = 1;
rate = accl(2,2) - accl(1,2);
for i = 1 : sz1
    % MJD time
    accl_mjd = accl(i,1);
    % Epoch's Day Comparison
    if fix(accl_mjd) - fix(mjd) == 0
        % Epoch's Time comparison
        accl_t = accl(i,2);
        if abs(t - accl_t) < 10^-7
            faccl = accl(i,3:5);
            facclx = faccl(1,1);
            faccly = faccl(1,2);
            facclz = faccl(1,3);
            break
        elseif t - accl_t < rate 
            [facclx] = interp_Lagrange(accl(i,:),accl(i,3),mjd,dpint);
            [faccly] = interp_Lagrange(accl(i,:),accl(i,4),mjd,dpint);
            [facclz] = interp_Lagrange(accl(i,:),accl(i,5),mjd,dpint);
            break
        end
    end
end
   