function [rgith,vgith] = dexxxeph_bcs2gcs(zbith,zbEMB,zgM, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dexxxeph_bcs2gcs : State vector transformation from Solar System 
%                    Barycentric coordinates to Geocentric coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the solar system bodies state vector with respect to
%  Geocenter. A simple transformation from Barycentric to Geocentric vector
%  is applied. The transformation is implemented through the Earth-Moon
%  barycentric vector.
%  
%
% Input arguments
% - zbith      : Barycentric state vector of the ith solar system item
% - zbEMB      : Barycentric state vector of the Earth-Moon barycenter
% - zgM        : Geocentric state vector of the Moon
% - HDfilename : DExxx header filename
%
%   >>       z = [x y z Vx Vy Vz]
%
% Output arguments:
% - zgith  : Geocentric state vector of the ith solar system item
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding Earth-Moon masses ratio according to the DExxx constants
%[GMconstant,AU,emrat,deformat,deperiod] = dexxxeph_readhd(HDfilename);
% Forces model settings matrix
planets_struct = orbit_model_struct.planets;

% GMconstant = planets_struct.DE_GMconstant;
% AU = planets_struct.DE_AU;
emrat = planets_struct.DE_EMRAT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Barycentric Position and velocity vetors
rbith = zbith(1,1:3)';
vbith = zbith(1,4:6)';
rbEMB = zbEMB(1,1:3)';
vbEMB = zbEMB(1,4:6)';
rgM = zgM(1,1:3)';
vgM = zgM(1,4:6)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geocentric Position and velocity vetors
rgith = rbith - (rbEMB - (1 / (1+emrat)) * rgM);
vgith = vbith - (vbEMB - (1 / (1+emrat)) * vgM);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
