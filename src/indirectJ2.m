function [ax,ay,az] = indirectJ2(C20,Re,GMmoon,rMoon,GMsun,rSun)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% indirectJ2 : Indirect J2 effect of Sun and Moon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  In case of Sun and Moon the interaction between the point-mass
%  gravitational attraction and the Earth's oblateness is taken account.
%  This is known as the so-called indirect J2 effect.
%
% Input arguments
% - C20     : full normalized second zonal harmonic coefficient
% - GMearth : Earth gravity constant  (m^3/sec^2)
% - Re      : radius  (meters)
% - GMmoon  : Moon gravity constant  (m^3/sec^2)
% - rMoon   : Moon body-fixed geocentric position vector (meters)
% - GMsun   : Sun gravity constant  (m^3/sec^2)
% - rSun    : Sun body-fixed geocentric position vector (meters)
%
% Output arguments:
% - ax,ay,az  :  Acceleration in Body-fixed system (m/s/s)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Celestial bodies Earth-centered coordinates and distance from Earth 
% Moon
xMoon = rMoon(1,1);
yMoon = rMoon(2,1);
zMoon = rMoon(3,1);
l_Moon = sqrt(rMoon(1,1)^2 + rMoon(2,1)^2 + rMoon(3,1)^2);
% Sun
xSun = rSun(1,1);
ySun = rSun(2,1);
zSun = rSun(3,1);
l_Sun = sqrt(rSun(1,1)^2 + rSun(2,1)^2 + rSun(3,1)^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Moon terms
termMoon = (GMmoon / l_Moon^3) * (Re / l_Moon)^2;
termMoon_x = termMoon * (5*(zMoon/l_Moon)^2 - 1) * xMoon;
termMoon_y = termMoon * (5*(zMoon/l_Moon)^2 - 1) * yMoon;
termMoon_z = termMoon * (5*(zMoon/l_Moon)^2 - 3) * zMoon;
% Sun terms
termSun = (GMsun / l_Sun^3) * (Re / l_Sun)^2;
termSun_x = termSun * (5*(zSun/l_Sun)^2 - 1) * xSun;
termSun_y = termSun * (5*(zSun/l_Sun)^2 - 1) * ySun;
termSun_z = termSun * (5*(zSun/l_Sun)^2 - 3) * zSun;

% Indirect J2 effect
indirectJ2_coef = -((3*sqrt(5))/2) * C20;

% Acceleration cartesian components
ax = indirectJ2_coef * (termMoon_x + termSun_x);
ay = indirectJ2_coef * (termMoon_y + termSun_y);
az = indirectJ2_coef * (termMoon_z + termSun_z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%