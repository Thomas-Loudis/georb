function [dstn,rms_stn,dz,rms_z,delta_kepler,rms_kepler,delta_Vstn rms_Vstn] = orbital_pert(z_ref,z_det,GM)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit comparison and peturbations transformation in the orbital frame
% 
% Purpose:
%  Comparison between reference and determined orbits (or different orbit
%  types) and peturbations transformation in the orbital frame (Gaussian).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remarks:
%  Transformation from inertial to orbital reference frame is based on the
%  rigorous mathematical scheme presented by Casotto (1993).
%
%  Reference (first) Orbit is used as the reference orbit in the final
%  equations of orbit peturbations.
%
%  Time is refered in Terrestrial Time (TT).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - z_ref: First (reference) orbit in celestial referene system GCRS
%          z_ref = [t r' v']
%          (e.g. PSO, Kimnematic Orbit)    
% - z_det: Second (determined) orbit in celestial referene system GCRS
%          z_det = [t r' v']
%          (e.g. Dynamic Orbit, Computed Orbit)
% - GM:    Earth gravity constant  (m^3/sec^2)
% 
% Output arguments:
% - dstn:    Orbit peturbations (or differences) in the orbital frame for
%            the three components, radial, along-track and cross-track
% - rms_stn: Orbit peturbations RMS in the orbital reference frame
%            rms_stn = [RMS_radial_i RMS_along_i RMS_cross_i]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                           May 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
%   06/05/2010  Upgrade for minimizing computation time (break). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if GM == 0
    GM = 3.986004415 * 10^14;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orbit differences are computed at common epochs of the two orbits
[j k] = size(z_det);
[m n] = size(z_ref);
epoch = 0;
wo = 1;
for i = 1 : j
    % Determined orbit
    t_det = z_det(i,1);
    r_det = [z_det(i,2) z_det(i,3) z_det(i,4)];
    v_det = [z_det(i,5) z_det(i,6) z_det(i,7)];
    % Orbit differences at common epochs
    for w = wo : m
        t_ref = z_ref(w,1);
        if abs(t_det - t_ref) < 10^(-8)
            wo = w + 1;
            % Common epoch of reference and determined orbits
            epoch = epoch + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Orbit differences in ICRS at the equivalent epoch
            r_ref = [z_ref(w,2) z_ref(w,3) z_ref(w,4)];
            v_ref = [z_ref(w,5) z_ref(w,6) z_ref(w,7)];         
            dt = t_det - t_ref;
            dr = r_det - r_ref;
            dv = v_det - v_ref;
            dz(epoch,:) = [dt dr dv];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            % Keplerian elements of Reference orbit z_ref (e.g. Kinematic)
            [a_k,e_k,incl_k,Omega_k,omega_k,f_k,M_k,E_k,u_k] = kepler(r_ref',v_ref',GM);
            kepler_ref = [a_k e_k incl_k Omega_k omega_k M_k];            
            %  conversion from degrees to radians
            incl_k = incl_k * (pi/180);
            Omega_k = Omega_k * (pi/180);
            omega_k = omega_k * (pi/180);
            f_k = f_k * (pi/180);
            M_k = M_k * (pi/180);
            E_k = E_k * (pi/180);
            u_k = u_k * (pi/180);
            % Keplerian elements of determined orbit z_det (e.g. Dynamic)
            [a_d,e_d,incl_d,Omega_d,omega_d,f_d,M_d,E_d,u_d] = kepler(r_det',v_det',GM);
            kepler_det = [a_d e_d incl_d Omega_d omega_d M_d];
            %  conversion from degrees to radians
            incl_d = incl_d * (pi/180) ;
            Omega_d = Omega_d * (pi/180);
            omega_d = omega_d * (pi/180);
            f_d = f_d * (pi/180);
            M_d = M_d * (pi/180);
            E_d = E_d * (pi/180);
            u_d = u_d * (pi/180);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Orbit differences of Keplerian elements (in degrees)
            delta_kepler(epoch,:) = kepler_det - kepler_ref;
            clear kepler_det kepler_ref
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Orbit differences in Orbital frame
            %   Components: Radial, Along-track, Cross-track
            %   Xsat = [ s t n ]'
            orbital_approach = 2;
            if (orbital_approach == 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st approach: Orbital Peturbations analysis             
            % Reference Orbit: First orbit z_ref (e.g. PSO,Kinematic orbit)
            rad_k = a_k * (1 - e_k * cos(E_k));
            n_ref = sqrt(GM / a_k^3);
            % Position perturbations
            d_radial(epoch,1) = (rad_k/a_k)*(a_d-a_k) - a_k*(e_d-e_k)*cos(f_k) + ((a_k*e_k)/sqrt(1-e_k^2))*(M_d-M_k)*sin(f_k);
            d_along(epoch,1) = a_k*(1+(1/(1-e_k^2))*(rad_k/a_k))*(e_d-e_k)*sin(f_k) + (a_k^2*sqrt(1-e_k^2)*(M_d-M_k))/rad_k + rad_k*((omega_d-omega_k)+(Omega_d-Omega_k)*cos(incl_k));
            d_cross(epoch,1)  = rad_k * ( (incl_d-incl_k)*sin(u_k) - (Omega_k-Omega_d)*sin(incl_k)*cos(u_k) );
            % Velocity perturbations
            d_Vradial(epoch,1) = -(n_ref * a_k * sin(f_k) / sqrt(1-e_k^2)) * (e_k*(a_d-a_k)/(2*a_k) + a_k*(e_d-e_k)/rad_k) - (n_ref*a_k^3/rad_k^2)*(M_d-M_k) - (n_ref*a_k^2*sqrt(1-e_k^2)/rad_k) * ((omega_d-omega_k)+(Omega_d-Omega_k)*cos(incl_k));
            d_Valong(epoch,1) = (n_ref*a_k / sqrt(1-e_k^2)) * ( (e_k+cos(f_k))*(e_d-e_k)/(1-e_k^2) - (1-e_k^2)*(a_d-a_k)/(2*rad_k) + e_k*((omega_d-omega_k)+(Omega_d-Omega_k)*cos(incl_k))*sin(f_k) );
            d_Vcross(epoch,1) = (n_ref*a_k / sqrt(1-e_k^2)) * ( (cos(u_k)+e_k*cos(omega_k))*(incl_d-incl_k) + sin(incl_k)*(sin(u_k)+e_k*sin(omega_k))*(Omega_d-Omega_k) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd apporach: Orbital frame transformation
            % Orbtial Frame transformatino matrix
            [Rrtn,er,et,en] = orbital_transf(r_ref',v_ref');

            delta_orbital_r = Rrtn * dr';
            d_radial(epoch,1) = delta_orbital_r(1);
            d_along(epoch,1)  = delta_orbital_r(2);
            d_cross(epoch,1)  = delta_orbital_r(3);
            
            delta_orbital_v = Rrtn * dv';
            d_Vradial(epoch,1) = delta_orbital_v(1);
            d_Valong(epoch,1)  = delta_orbital_v(2);
            d_Vcross(epoch,1)  = delta_orbital_v(3);            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
            tcommon(epoch,1) = t_ref;
            break
        end        
    end    
end
clear j k m n i w 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit differences/peturbations matrices
dstn = [tcommon d_radial d_along d_cross];
dz = [tcommon dz(:,2:7)];
delta_kepler = [tcommon delta_kepler];
delta_Vstn = [tcommon d_Vradial d_Valong d_Vcross];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbital Frame
[sz5 sz6] = size(dstn);
rms_stn = zeros(1,3);
rms_Vstn = zeros(1,3);
j = 1;
for i = 2 : sz6
    rms_stn(1,j) = rms(dstn(:,i));
    rms_Vstn(1,j) = rms(delta_Vstn(:,i));
    j = j + 1;
end
clear sz5 sz6 j i
% Celestial Frame and Kepler elements
[sz5 sz6] = size(dz);
rms_z = zeros(1,6);
rms_kepler = zeros(1,6);
j = 1;
for i = 2 : sz6
    rms_z(1,j) = rms(dz(:,i));
    rms_kepler(1,j) = rms(delta_kepler(:,i));
    j = j + 1;
end
clear sz5 sz6 j i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
