function [a_relativistic] = relativistic_effects(Zsat,Zearth,GMearth,GMsun,cslight,ppn_beta,ppn_gama,effect_name)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relativistic corrections to satellite's acceleration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the relativistic corrections to the acceleration of an
%  Earth's artificial satellite.
%
%  The formula used here is in agreement with IERS Conventions 2010 which
%  is based to formalisms of Brumberg and Kopeikin, Damour, Soffel and Xu 
%  and the IAU Resolutions B1.3 and B1.4 (2000).
%
%  These formalisms allow a full post-Newtonian framework 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - Zsat        : Satellite's state vector with respect to the Earth
% - Zearth      : Earth's state vector with respect to the Sun 
%               >> z = [x y z Vx Vy Vz]'
% - GMearth     : Gravitational coefficient of Earth 
% - GMsun       : Gravitational coefficient of Sun
% - cslight     : Speed of light
% - ppn_beta    : PPN (parameterized post-Newtonian) parameters
% - ppn_gama    : PPN equal to 1 in General Relativity
% - effect_name : Relativistic effect name to be computed:
%                 'Schwarzschild' or 
%                 'Lense-Thirring' or 
%                 'deSitter'  referring to the geodetic effect of de Sitter precession
%
% Output arguments:
% - a_relativistic  : Relativistic effect corrections to acceleration vector
%                     a_relativistic = [ax; ay; az]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         June 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
%  9/12/2022    Dr. Thomas Loudis Papanikolaou,
%               Minor Code upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite vectors
rsat = [Zsat(1,1) Zsat(2,1) Zsat(3,1)]';
vsat = [Zsat(4,1) Zsat(5,1) Zsat(6,1)]';
% Earth-Satellite distance
lsat = sqrt(rsat(1,1)^2 + rsat(2,1)^2 + rsat(3,1)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth vectors
rearth = [Zearth(1,1) Zearth(2,1) Zearth(3,1)]';
vearth = [Zearth(4,1) Zearth(5,1) Zearth(6,1)]';

% Earth's angular momentum per unit mass
Jearth_L = productcross(rearth,vearth);
GM_Earth = 3.986004418 * 10^14;
G = 6.67428 * 10^-11;
M_Earth = GM_Earth / G;
Jearth = Jearth_L / M_Earth;
% Earth's angular momentum per unit mass (m^2/sec) as per IERS Conv. 2010
Jearth = [0; 0; 9.8*10^8];              

% Sun-Earth distance
learth = sqrt(rearth(1,1)^2 + rearth(2,1)^2 + rearth(3,1)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relativistic effects terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwarzschild effect
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = strcmp(effect_name,'Schwarzschild');
%if relativ(1,1) == 1
if test == 1
% Schwarzschild terms
a_Schwarzschild = (GMearth / (cslight^2 * lsat^3)) * ...
                  ( ...
                  ( 2*(ppn_beta+ppn_gama) * (GMearth/lsat) - ppn_gama * productdot(vsat,vsat) ) * rsat ...
                  + 2 * (1 + ppn_gama) * productdot(rsat,vsat) * vsat ...
                  );
else
    a_Schwarzschild = [0; 0; 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lense-Thirring precession (frame-dragging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = strcmp(effect_name,'Lense-Thirring');
if test == 1
%if relativ(1,2) == 1
a_LenseThirring = (1 + ppn_gama) * (GMearth / (cslight^2 * lsat^3)) * ...
                  ( (3/lsat^2) * productcross(rsat,vsat) * productdot(rsat,Jearth) + productcross(vsat,Jearth) );              
else
    a_LenseThirring = [0; 0; 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geodetic effect or de Sitter precession, geodesic precession
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test = strcmp(effect_name,'deSitter');
if test == 1
%if relativ(1,3) == 1
vec1 = - ( GMsun / (cslight^2 * learth^3) ) * rearth;
cp1 = productcross(vearth,vec1);
cp2 = productcross(cp1,vsat);
a_deSitter = (1 + 2 * ppn_gama) * cp2;
%clear vec1 cp1 cp2
else
    a_deSitter = [0; 0; 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Vectors magnitude
% a_Schwarzschild_magn = sqrt(a_Schwarzschild(1,1)^2 + a_Schwarzschild(2,1)^2 + a_Schwarzschild(3,1)^2)
% a_LenseThirring_magn = sqrt(a_LenseThirring(1,1)^2 + a_LenseThirring(2,1)^2 + a_LenseThirring(3,1)^2)
% a_deSitter_magn      = sqrt(a_deSitter(1,1)^2 + a_deSitter(2,1)^2 + a_deSitter(3,1)^2)

% Overall vector
a_relativistic = a_Schwarzschild + a_LenseThirring + a_deSitter;
