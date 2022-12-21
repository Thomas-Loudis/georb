function [partials_2nd_rpl, partials_2nd_xyz, partials_1st_rpl, partials_1st_xyz] = harmonics_partials_2nd(r,n_max,m_max,GM,ae,Cnm,Snm)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd order Partial derivatives of the gravitational potential 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Computation of the 2nd order partial derivatives of the gravitational
%   potential w.r.t. spherical and cartestian coordinates
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The partial derivatives are noted as:
%   dfx / dy = Uxy,     dfx / dx = Uxx
% and they are related 
%   Uxy = Uyx, Uxz = Uzx, Uyz = Uzy and
% Laplace's equation or Laplacian:
%   Uxx + Uyy + Uzz = 0 
%
% Matrix of second partial derivatives w.r.t. Cartesian coordinates:
%   U = [ Uxx   Uxy   Uxz 
%         Uxy   Uyy   Uyz
%         Uxz   Uyz   Uzz ]  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - r:                  Position vector (m)
% - GM:                 Earth gravity constant  (m^3/sec^2)
% - ae:                 radius  (meters)
% - Cnm, Snm:           normalized spherical harmonics coefficients
% - n_max:              Truncation Degree of harmoncis series expansion 
% - m_max:              Truncation Order of harmoncis series expansion 
%
% Output arguments:
% - partials_1st_rpl:   1st order partials w.r.t. spherical coordinates
%                       radius, latitude, longitude
% - partials_1st_xyz:   1st order partials w.r.t. Cartesian coordinates xyz
% - partials_2nd_rpl:   2nd order partials w.r.t. spherical coordinates
%                       radius, latitude, longitude
% - partials_2nd_xyz:   2nd order partials w.r.t. Cartesian coordinates xyz
%                       radius, latitude, longitude
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         June 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 09/12/2022  Thomas Loudis Papanikolaou
%             MInor code upgrade and rename as function harmonics_partials_2nd.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st order of partials of gravitational potential
[partials_1st_rpl, partials_1st_xyz] = harmonics_partials_1st(r,n_max,m_max,GM,ae,Cnm,Snm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dV_r     = partials_1st_rpl(1,1);
dV_phi   = partials_1st_rpl(2,1);
dV_lamda = partials_1st_rpl(3,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of spherical coordinates (in radians)
[lamda,phi,l] = lamda_phi(r);
rdist = l;
% Normalized associated Legendre functions
[Pnm_norm] = Legendre_functions(phi,n_max);

% First-order derivatives of normalized associated Legendre functions
[dPnm_norm] = Legendre1ord(phi,n_max) ;

% Second-order derivatives of the Normalized Associated Legendre functions
[d2Pnm_norm] = Legendre2ord(phi,n_max) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geopotential 2nd-order partial derivatives in radius,theta,lamda components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vrr, VrTheta (Vrt), VrLamda(VrL), V_ThetaTheta (Vtt), V_ThetaLamda (VtL)
% V_LamdaLamda (VLL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vrr = 0; Vrt = 0; VrL = 0; Vtt = 0; VtL = 0; VLL = 0;
for n = 0 : n_max
    for m = 0 : n
        Vrr_f = (GM / ae^3) * (n+1)*(n+2) * (ae / l)^(n+3) ;
        Vrr = Vrr + Vrr_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * Pnm_norm(n+1,m+1) ;
        
        Vrt_f = - (GM / ae^2) * (n+1) * (ae / l)^(n+2);  % Eq.6.24 Revised 2012
        % Abramowitz and Stegun 1972 (first-order derivative of Pnm)
        %dPnm = n*tan(phi) * Pnm_norm(n+1,m+1) - (n+m) * (1/cos(phi)) * Pnm_norm(n-1+1,m+1)  ;
        Vrt = Vrt + Vrt_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * dPnm_norm(n+1,m+1) ;
        %Pnm_dv1 = dPnm_norm(n+1,m+1)
        
        VrL_f = (GM / ae^2) * (n+1) * (ae / l)^(n+2) ;
        VrL = VrL + VrL_f * m * ( Cnm(n+1,m+1) * sin(m*lamda) - Snm(n+1,m+1) * cos(m*lamda) ) * Pnm_norm(n+1,m+1) ;
        
        Vtt_f = (GM / ae) * (ae / l)^(n+1) ;
        Vtt = Vtt + Vtt_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * d2Pnm_norm(n+1,m+1);
        % Pnm_dv2 = d2Pnm_norm(n+1,m+1)
        % % Abramowitz and Stegun 1972 (first-order derivative of Pnm)
        % % dPnm = n*tan(phi) * Pnm_norm(n+1,m+1) - (n+m) * (1/cos(phi)) * Pnm_norm(n-1+1,m+1)  ;
        % % Hobson 1931 (second-order derivative of Pnm)
        % dPnm = dPnm_norm(n+1,m+1);
        % dPnm2 = - tan(phi) * dPnm - ( n*(n+1) - m^2 * (1/cos(phi))^2 ) * Pnm_norm(n+1,m+1) ;
        % Vtt = Vtt + Vtt_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * dPnm2;  %d2Pnm_norm(n+1,m+1) ;

        VtL_f = (GM / ae) * (ae / l)^(n+1) ;
        % Corrected : May 2012
        VtL = VtL + VtL_f * m * ( - Cnm(n+1,m+1) * sin(m*lamda) + Snm(n+1,m+1) * cos(m*lamda) ) * dPnm_norm(n+1,m+1) ;

        VLL_f = (GM / ae) * (ae / l)^(n+1) ;
        VLL = VLL - VLL_f * m^2 * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * Pnm_norm(n+1,m+1) ;
    end
end

% Geopotential 2nd-order partial derivatives in local tangent frame [e_radius, e_theta, e_lamda]
Vrtl = [ Vrr Vrt VrL
         Vrt Vtt VtL
         VrL VtL VLL ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives of (r,phi,lamda) with respect to (x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDVrx = [
%       cos(phi)*cos(lamda)            cos(phi)*sin(lamda)          sin(phi)
%   (-1/l)*sin(phi)*cos(lamda)     (-1/l)*sin(phi)*sin(lamda)    (1/l)*cos(phi)
% ( -1/(l*cos(phi)) )*sin(lamda)  ( 1/(l*cos(phi)) )*cos(lamda)        0
% ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives of (r,theta,lamda) with respect to (x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDVrx = [
%       sin(theta)*cos(lamda)            sin(theta)*sin(lamda)           cos(theta)
%   ( 1/l)*cos(theta)*cos(lamda)     ( 1/l)*cos(theta)*sin(lamda)    (-1/l)*sin(theta)
% ( -1/(l*sin(theta)) )*sin(lamda)  ( 1/(l*sin(theta)) )*cos(lamda)         0
% ];
% Replacement of "theta" with "phi"
PDVrx = [
      cos(phi)*cos(lamda)            cos(phi)*sin(lamda)          sin(phi)
  ( 1/l)*sin(phi)*cos(lamda)     ( 1/l)*sin(phi)*sin(lamda)    (-1/l)*cos(phi)
( -1/(l*cos(phi)) )*sin(lamda)  ( 1/(l*cos(phi)) )*cos(lamda)        0
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Theta angle
theta = pi/2 - phi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives of PDVrx with respect to (r,theta,lamda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdv2_r = [
              0                                      0                                 0
(-1/rdist^2)*(cos(theta)*cos(lamda))  (-1/rdist^2)*(cos(theta)*sin(lamda))   (1/rdist^2)*sin(theta)
(1/rdist^2)*(sin(lamda)/sin(theta))   (-1/rdist^2)*(cos(lamda)/sin(theta))             0
];    

pdv2_theta = [
             cos(theta)*cos(lamda)                             cos(theta)*sin(lamda)                  -sin(theta)
       (-1/rdist)*(sin(theta)*cos(lamda))              (-1/rdist)*(sin(theta)*sin(lamda))         (-1/rdist)*cos(theta)
(1/rdist)*(sin(lamda)*cos(theta)/sin(theta)^2)  (-1/rdist)*(cos(lamda)*cos(theta)/sin(theta)^2)             0
];

pdv2_lamda = [
       -sin(theta)*sin(lamda)               sin(theta)*cos(lamda)         0
(-1/rdist)*(cos(theta)*sin(lamda))   (1/rdist)*(cos(theta)*cos(lamda))    0
(-1/rdist)*(cos(lamda)/sin(theta))  (-1/rdist)*(sin(lamda)/sin(theta))    0
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives of (Fx,Fy,Fz) with respect to (r,theta,lamda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pdv(fxyz / r)
pdvfx_r = Vrr * PDVrx(1,1) + Vrt * PDVrx(2,1) + VrL * PDVrx(3,1) + ...
          dV_r * pdv2_r(1,1) + dV_phi * pdv2_r(2,1) + dV_lamda * pdv2_r(3,1);
      
pdvfy_r = Vrr * PDVrx(1,2) + Vrt * PDVrx(2,2) + VrL * PDVrx(3,2) + ...
          dV_r * pdv2_r(1,2) + dV_phi * pdv2_r(2,2) + dV_lamda * pdv2_r(3,2);
      
pdvfz_r = Vrr * PDVrx(1,3) + Vrt * PDVrx(2,3) + VrL * PDVrx(3,3) + ...
          dV_r * pdv2_r(1,3) + dV_phi * pdv2_r(2,3) + dV_lamda * pdv2_r(3,3);

% pdv(fxyz / theta)
pdvfx_theta = Vrt * PDVrx(1,1) + Vtt * PDVrx(2,1) + VtL * PDVrx(3,1) + ...
              dV_r * pdv2_theta(1,1) + dV_phi * pdv2_theta(2,1) + dV_lamda * pdv2_theta(3,1);
      
pdvfy_theta = Vrt * PDVrx(1,2) + Vtt * PDVrx(2,2) + VtL * PDVrx(3,2) + ...
              dV_r * pdv2_theta(1,2) + dV_phi * pdv2_theta(2,2) + dV_lamda * pdv2_theta(3,2);

pdvfz_theta = Vrt * PDVrx(1,3) + Vtt * PDVrx(2,3) + VtL * PDVrx(3,3) + ...
              dV_r * pdv2_theta(1,3) + dV_phi * pdv2_theta(2,3) + dV_lamda * pdv2_theta(3,3);

% pdv(fxyz / lamda)
pdvfx_lamda = VrL * PDVrx(1,1) + VtL * PDVrx(2,1) + VLL * PDVrx(3,1) + ...
              dV_r * pdv2_lamda(1,1) + dV_phi * pdv2_lamda(2,1) + dV_lamda * pdv2_lamda(3,1);

pdvfy_lamda = VrL * PDVrx(1,2) + VtL * PDVrx(2,2) + VLL * PDVrx(3,2) + ...
              dV_r * pdv2_lamda(1,2) + dV_phi * pdv2_lamda(2,2) + dV_lamda * pdv2_lamda(3,2);
          
pdvfz_lamda = VrL * PDVrx(1,3) + VtL * PDVrx(2,3) + VLL * PDVrx(3,3) + ...
              dV_r * pdv2_lamda(1,3) + dV_phi * pdv2_lamda(2,3) + dV_lamda * pdv2_lamda(3,3);
          
% matrix : pdvFxyz_rtl
pdvFxyz_rtl = [
pdvfx_r   pdvfx_theta   pdvfx_lamda
pdvfy_r   pdvfy_theta   pdvfy_lamda
pdvfz_r   pdvfz_theta   pdvfz_lamda
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geopotential 2nd-order partial derivatives in ITRS X,Y,Z components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Umatrix = PDVrx' * pdvFxyz_rtl';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

partials_2nd_rpl = Vrtl;
partials_2nd_xyz = Umatrix;
