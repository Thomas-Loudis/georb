function [orbc] = orbt2c(orbt,eopdat,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbt2c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Orbit transformation from the Terrestrial Reference Frame (ITRF) to the
%  Celestial Reference Frame (GCRF).  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - orbt:   Orbit in ITRF  (i.e. orbt(i,:) = [mjd rITRF' vITRF'] )
% - eopdat: Earth Orientation Parameters (EOP) data that are required for
%           the orbit arc length
% - dpint:  Number of data points (days) that are required for the EOP
%           interpolation to the computation epoch
%
% Output arguments:
% - orbc:   Orbit in GCRS  (i.e. orbt(i,:) = [mjd rGCRS' vGCRS'] )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   April 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 08/07/2012   Distinguishing for Velocity vector (incuding or not)
%              In case of position vector only (Kinematic orbit data)
%              transformation is applied only to position vector coords.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for Velocity vector
[d1 d2] = size(orbt);
if d2 < 7
    Vel_vector = 0;
elseif d2 >= 7
    Vel_vector = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit transformation from ITRS to GCRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz1 sz2] = size(orbt);
if Vel_vector == 1
    orbc = zeros(sz1,7);
elseif Vel_vector == 0
    orbc = zeros(sz1,4);
end
for it = 1 : sz1
    [EOP,dEOP] = trs2crs(orbt(it,1),eopdat,dpint);
    rGCRS = (EOP) * orbt(it,2:4)';
    if Vel_vector == 1
        vGCRS = (EOP) * orbt(it,5:7)' + dEOP * orbt(it,2:4)';
        orbc(it,:) = [orbt(it,1) rGCRS' vGCRS'];
    elseif Vel_vector == 0
        orbc(it,:) = [orbt(it,1) rGCRS'];
    end
end
end