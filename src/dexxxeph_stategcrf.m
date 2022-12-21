function [rg3rd_meter,vg3rd_meter] = dexxxeph_stategcrf(body_ith,jd,HDxxxfilename,DEcheby,DErecord)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dexxxeph_stategcrf :  Computation of Celestial body's state vector in GCRF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the state vector for Sun, Moon and Planets in units of
%  meters and m/sec in the Geocentric Celestial Reference Frame.
%  Coordinates computation is based on the JPL's DExxx series ephemerides
%  and the evaluation of Chebychev polynomials approximation.
%
% Input arguments
% - body_ith   : Number of the solar system body according to the DExxx
%                list of the 13 items
% - tjd        : Julian Date number of the required epoch
% - HDfilename : Header file name of the DExxx ephemeris
% - decbv      : DExxx Chebyshev coefficients for the required data record
% - derecord   : DExxx data record's number, start, end JD time
%   >>           derecord = [ datarecord_no JDstart JDend]'
%
% Output arguments:
% - r         : Position vector in Solar Barycentric System   (Km)
% - v         : Velocity vector in Solar Barycentric System   (Km/day)
% - chebycoef : Chebyshev coefficients that were required only
%
% >> r = [ X Y Z]';  v = [Vx Vy Vz]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Moon Geocentric coordinates (Km,Km/d)
[rgMoon,vgMoon,chebycoef] = dexxxeph_state(10,jd,HDxxxfilename,DEcheby,DErecord);
zgMoon = [rgMoon' vgMoon'];

% Earh-Moon-Barycenter coordinates in Solar System Barycentrer (Km,Km/d)
[rEMB,vEMB,chebycoef] = dexxxeph_state(3,jd,HDxxxfilename,DEcheby,DErecord);
zbEMB = [rEMB' vEMB'];

% Planetary Coordinates
if body_ith == 10
    rg3rd = rgMoon;
    vg3rd = vgMoon;
else
    % Celestial body's state vector
    [rb3rd,vb3rd,chebycoef] = dexxxeph_state(body_ith,jd,HDxxxfilename,DEcheby,DErecord);
    zb3rd = [rb3rd' vb3rd'];
    % Transformation : Barycentric to Geocentric Celestial Reference Frame
    [rg3rd,vg3rd] = dexxxeph_bcs2gcs(zb3rd,zbEMB,zgMoon,HDxxxfilename);
end

% State vector Units Conversion: Km and Km/day to m and m/sec
rg3rd_meter = rg3rd * 10^3;
vg3rd_meter = vg3rd * (10^3 / (24*3600));
