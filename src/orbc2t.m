function [orbt] = orbc2t(orbc,eopdat,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbc2t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Orbit transformation from the Celestial Reference Frame (GCRF) to the 
%  Terrestrial Reference Frame (ITRF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - orbc:   Orbit in GCRS  (i.e. orbc(i,:) = [mjd rGCRS' vGCRS'] )
% - eopdat: Earth Orientation Parameters (EOP) data that are required for
%           the orbit arc length
% - dpint:  Number of data points (days) that are required for the EOP
%           interpolation to the computation epoch
%
% Output arguments:
% - orbt:   Orbit in ITRS  (i.e. orbt(i,:) = [mjd rITRS' vITRS'] )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   April 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 13/05/2021, Dr. Thomas Papanikolaou
%             Test for the velocity vector of the input orbit array  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity vector
% (Kinematic Orbit Data may not include velocity)
[d1, d2] = size(orbc);
if d2 < 7
    Vel_vector = 0;
elseif d2 >= 7
    Vel_vector = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix preallocation
orbt = zeros(d1,d2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit transformation from GCRS to ITRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it = 1 : d1
    mjd_it = orbc(it,1);
    [EOP,dEOP,EOP_inv,dEOP_inv] = trs2crs(mjd_it,eopdat,dpint);
    rITRS = (EOP)' * orbc(it,2:4)';
    if Vel_vector == 0        
    orbt(it,:) = [mjd_it rITRS'];
    elseif Vel_vector == 1
    vITRS = (EOP)' * orbc(it,5:7)' + dEOP_inv * orbc(it,2:4)';
    orbt(it,:) = [mjd_it rITRS' vITRS'];
    end
end