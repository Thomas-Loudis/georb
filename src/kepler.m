function [a,e,i,Omega,omega,f,M,E,u] = kepler(r,v,GM)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keplerian elements from Position and Velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Conversion from state vector (z) components to the Keplerian elements
%  [a,e,i,Omega,omega,f,M,E,u]
%
% Input arguments:
% - t:  seconds since 0 hours is refered in TT (Terrestrial Time)
% - r:  position in Celestial Reference System (GCRS), m
% - v:  velocity in Celestial Reference System (GCRS), m/sec
% - GM: Earth gravity constant,  m^3/sec^2
%
% Output arguments:
%  - a:      semi-major axis  (m)
%  - e:      eccentricity
%  - i:      inclination
%  - Omega:  right ascension of the ascending node
%  - omega:  argument of perigee
%  - f:      true anomaly
%  - M:      Mean anomaly
%  - E:      eccentric anomaly
%  - u:      argument of latitude
%
%  The angular elements are computed in degrees.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                      November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% position vector r = [ X Y Z ]'
X = r(1,1);
Y = r(2,1);
Z = r(3,1);
% velocity vector v = [ Vx Vy Vz ]'
Vx = v(1,1);
Vy = v(2,1);
Vz = v(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h: angular momentum (areal velocity vector)
h = [ Y*Vz-Z*Vy
      Z*Vx-X*Vz
      X*Vy-Y*Vx ];

W = h/sqrt(h'*h);

Wx = W([1]);
Wy = W([2]);
Wz = W([3]);

% parameter p: semi-latus rectum
p = (h'*h)/GM;

% i: the inclination in radians
i = arctan( sqrt(Wx^2+Wy^2) , Wz );
% converse in degrees
i_deg = i * (180/pi);

% Omega: the right ascension of the ascending node in radians
Omega = arctan( Wx , -Wy );
% converse in degrees
Omega_deg = Omega * (180/pi);

% a: the semi-major axis
a = 1/( (2/sqrt(r'*r))-(v'*v)/ GM );

% n: the mean motion
n = sqrt( GM /a^3);

% e: the eccentricity
e = sqrt(1- (p/a) );

% E: the eccentric anomaly in radians
E = arctan( (r'*v)/(a^2*n) , 1-sqrt(r'*r)/a );
% converse in degrees
E_deg = E * (180/pi);

% M: the mean anomaly in radians
M = E - e*sin(E);
% converse in degrees
M_deg = M * (180/pi);

% u: the argument of latitude in radians
u = arctan( Z , -X*Wy+Y*Wx );
% converse in degrees
u_deg = u * (180/pi);

% f : the true anomaly in radians
f = arctan( sqrt(1-e^2)*sin(E) , cos(E)-e );
% converse in degrees
f_deg = f * (180/pi);

% omega : the argument of perigee in degrees
omega_deg = u_deg - f_deg;
if omega_deg < 0
    omega_deg = omega_deg + 360;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  The angular elements are computed in degrees.
i = i_deg;
Omega = Omega_deg;
omega = omega_deg;
f = f_deg;
E = E_deg;
M = M_deg;
u = u_deg;
