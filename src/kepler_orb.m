function [z,keplerian] = kepler_orb(zo,tmax,h,GM)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference Keplerian Orbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of state vector from "orbit design" Keplerian elements at
%  equidistant points.
%
%  Equidistant epochs are defined according to the selected integration
%  step "h" and the orbit arc length.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - zo: Initial epoch and state vector in Celestial Reference System GCRS
%       zo = [to ro' vo']
%       to: initial epoch in seconds of TT (Terrestrial Time) since 0 (h)
%       ro: initial position in GCRS (m)
%       vo: initial velocity in GCRS (m/sec)
%
% Output arguments:
% - z: State vector
%      zi = [ti Xi Yi Zi Vxi Vyi Vzi], for individual epoch "ti"
%   t: epoch since 0h in seconds 
%   r: position in Celestial Reference System (GCRS), m
%   v: velocity in Celestial Reference System (GCRS), m/sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         May 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial epoch
to = zo(1,1);
% Initial position and velocity vectors in the GCRS
ro = [zo(1,2); zo(1,3); zo(1,4)];
vo = [zo(1,5); zo(1,6); zo(1,7)];
% Store initial epoch State Vector
z(1,:) = zo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kepler elements at initial epoch (to)
[a,e,incl,Omega,omega,fo,Mo,Eo,uo] = kepler(ro,vo,GM);
% Store Keplerian elements
keplerian(1,:) = [a e incl Omega omega Mo];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  converse in radians
Mo_rad = Mo * (pi/180);

% mean motion n
n = sqrt( GM /a^3);

indx = 2;
for t = to+h : h : tmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mean anomaly
    M_rad = Mo_rad + n * (t - to);
    % converse in degrees
    Mdeg = M_rad * (180/pi);
    nrev = fix(Mdeg /  360);
    if nrev >=1
        M = Mdeg - nrev * 360;
    else
        M = Mdeg;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Eccentric anomaly
    [E] = kepler_eq(M,e);
    % Keplerian elements per epoch matrix
    kepl = [a e incl Omega omega E];
    % Store Keplerian elements
    keplerian(indx,:) = [a e incl Omega omega M];
    % State vector
    [z_kepl] = kepler_z(kepl);
    % Store next epoch State Vector
    z(indx,:) = [t z_kepl];
    indx = indx + 1;
end