function [z] = kepler_z(kepler)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position and Velocity from Keplerian elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Conversion from the 6 Keplerian elements (a,e,i,Omega,omega,E) to the
%  state vector (z).
%
%  Input/Output arguments can be either values or matrices.
%
% Input arguments:
%  - kepler:  Keplerian elements
%             kepler = [a e i Omega omega E]
%    a:       semi-major axis                        (m)
%    e:       eccentricity
%    i:       inclination
%    Omega:   right ascension of the ascending node
%    omega:   argument of perigee
%    E:       eccentric anomaly
% -  GM:      Earth gravity constant                 (m^3/sec^2)
%
%  The angular elements of the input arguments are in degrees.
%
% Output arguments:
% - z:  State vector
%       zi = [Xi Yi Zi Vxi Vyi Vzi], for individual epoch "i"
%   r:  position in Celestial Reference System (GCRS), m
%   v:  velocity in Celestial Reference System (GCRS), m/sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                     November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 8/11/2022    Dr. Thomas Papanikolaou
%              Code minor upgarde 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
a = kepler(1,1);
e = kepler(1,2);
% conversion of i from degrees to radians
i = deg2rad(kepler(1,3));

% Conversion of Omega from degrees to radians
Omega = deg2rad(kepler(1,4));

% conversion of omega from degrees to radians
omega = deg2rad(kepler(1,5));

% conversion of E from degrees to radians
E = deg2rad(kepler(1,6));

% GM: Earth gravitational constant in m^3/sec^2 by standards of IERS Conventions 2003
GM = 3986004.418*10^8;   

% n : mean motion
n = sqrt(GM/a^3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = [cos(omega)*cos(Omega)-sin(omega)*cos(i)*sin(Omega)
     cos(omega)*sin(Omega)+sin(omega)*cos(i)*cos(Omega)
     sin(omega)*sin(i)];

T = [-sin(omega)*cos(Omega)-cos(omega)*cos(i)*sin(Omega)
     -sin(omega)*sin(Omega)+cos(omega)*cos(i)*cos(Omega)
      cos(omega)*sin(i)];

N = [sin(i)*sin(Omega) 
    -sin(i)*cos(Omega)
     cos(i)];

% R : rotation matrix
R = [S T N];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position vector on the orbital plane
X = [a*(cos(E)-e); a*sqrt(1-e^2)*sin(E); 0;];

% Velocity vector on the orbital plane
dX = ((n*a)/(1-e*cos(E)))*[-sin(E); sqrt(1-e^2)*cos(E); 0;];

% r : position vector  
r = R*X;

% v : velocity vector
v = R*dX;

% z : state vector
z_vec = [r; v;];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = z_vec';
