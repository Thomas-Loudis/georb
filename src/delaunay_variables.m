function [F1,F2,F3,F4,F5] = delaunay_variables(MJD)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: delaunay_variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
%  Computation of the Delaunay variables according to the IERS Conventions
%  2010.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - mjd             : MJD including the fraction of the day
%
% Output arguments:
% - F1,F2,F3,F4,F5  : Delaunay variables (l,l',F,D,Omega) in radians
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                    June 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficient for Conversion from arcsec to radians
arcsec2rad = pi / (180 * 3600);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter t
JD_TT = MJD + 2400000.5;
taph = ( JD_TT - 2451545.0 ) / 36525;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delaunay variables for Sun and Moon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F1 = l = Mean anonaly of the Moon
F1 = 134.96340251 * (pi/180) + 1717915923.2178 * arcsec2rad * taph + 31.8792 * arcsec2rad * taph^2 + 0.051635 * arcsec2rad * taph^3 - 0.00024470 * arcsec2rad * taph^4;
% F2 = l' = Mean anonaly of the Sun
F2 = 357.52910918 * (pi/180) + 129596581.0481 * arcsec2rad * taph - 0.5532 * arcsec2rad * taph^2 + 0.000136 * arcsec2rad * taph^3 - 0.00001149 * arcsec2rad * taph^4;
% F3 = F = L - Omega
F3 = 93.27209062 * (pi/180) + 1739527262.8478 * arcsec2rad * taph - 12.7512 * arcsec2rad * taph^2 - 0.001037 * arcsec2rad * taph^3 + 0.00000417 * arcsec2rad * taph^4;
% F4 = D = Mean Elongation of the Moon from the Sun
F4 = 297.85019547 * (pi/180) + 1602961601.2090 * arcsec2rad * taph - 6.3706 * arcsec2rad * taph^2 + 0.006593 * arcsec2rad * taph^3 - 0.00003169 * arcsec2rad * taph^4;
% F5 = Omega = Mean Longitude of the Ascending Node of the Moon
F5 = 125.04455501 * (pi/180) - 6962890.5431 * arcsec2rad * taph + 7.4722 * arcsec2rad * taph^2 + 0.007702 * arcsec2rad * taph^3 - 0.00005939 * arcsec2rad * taph^4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
