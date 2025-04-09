function [partials_param, partials_c, partials_s] = harmonics_partials_coef(r,n_max,n_min,GM,ae,gravity_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial derivatives of the gravitational potential w.r.t. gravity field parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the partial derivatives of the gravitational potential
%  w.r.t. spherical harmonics coefficients
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
% - partials_c:     partials w.r.t. harmonics coefficients Cnm 
% - partials_s:     partials w.r.t. harmonics coefficients Snm 
% - partials_p:     OVerall matrix of partials w.r.t. harmonics coefficients Cnm and Snm 
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                  22 March 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field model structure matrix
gfm_struct = gravity_struct;
C_degree_order = gfm_struct.C_degree_order_estim;
S_degree_order = gfm_struct.S_degree_order_estim;
% Gravity Field parameters :: C,S coefficients' degree and order
[Nparam_C, k] = size(C_degree_order);
[Nparam_S, k] = size(S_degree_order);
degree_max = C_degree_order(Nparam_C,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of spherical coordinates in radians
[lamda,phi,l] = lamda_phi(r);      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation of normalized associated Legendre functions
[Pnm_norm] = Legendre_functions(phi,n_max);
% First-order derivatives of normalized associated Legendre functions
[dPnm_norm] = Legendre1ord(phi,n_max) ;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices Pre-allocation
% Nrange_coef = (n_max - n_min) + 1;
% Ncoef = 0;
% for n = n_min : n_max
%     for m = 0 : n
%         Ncoef = Ncoef + 1;
%     end
% end

Nparam_Cnm = Nparam_C;
Nparam_Snm = Nparam_S;
partials_c = zeros(3 , Nparam_Cnm);
partials_s = zeros(3 , Nparam_Snm);
partials_param = zeros(3 , Nparam_Cnm+Nparam_Snm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if 1 < 0
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_coef = 0;
for n = n_min : n_max
    for m = 0 : n
        
%         dF_r_Cnm     = -(GM / l^2) * (n+1) * (ae/l)^(n+2) * Pnm_norm(n+1,m+1) * cos(m*lamda);  
%         dF_theta_Cnm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * cos(m*lamda) ; 
%         dF_lamda_Cnm =  (GM / l)   * (ae/l)^n * m * ( -sin(m*lamda) ) * Pnm_norm(n+1,m+1);  
% 
%         dF_r_Snm     = -(GM / l^2) * (n+1) * (ae/l)^(n+2) * Pnm_norm(n+1,m+1) * sin(m*lamda);
%         dF_theta_Snm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * sin(m*lamda) ;
%         dF_lamda_Snm =  (GM / l)   * (ae/l)^n * m * cos(m*lamda) * Pnm_norm(n+1,m+1) ; 

        % dV_r = dV_r + (GM/l^2) * ( -(n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) * (Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda)) );
        % dV_phi = dV_phi     + (GM / l) * ((ae/l)^n) * dPnm_norm(n+1,m+1) * (Cnm(n+1,m+1)*cos(m*lamda)+Snm(n+1,m+1)*sin(m*lamda));
        % dV_lamda = dV_lamda + (GM / l) * m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (Snm(n+1,m+1)*cos(m*lamda)-Cnm(n+1,m+1)*sin(m*lamda));
        
        dF_r_Cnm     = -(GM / l^2) * (n+1) * (ae/l)^n * Pnm_norm(n+1,m+1) * cos(m*lamda); 
        dF_theta_Cnm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * cos(m*lamda) ; 
        dF_lamda_Cnm =  (GM / l)   * (ae/l)^n * m * ( -sin(m*lamda) ) * Pnm_norm(n+1,m+1);  

        dF_r_Snm     = -(GM / l^2) * (n+1) * (ae/l)^n * Pnm_norm(n+1,m+1) * sin(m*lamda);
        dF_theta_Snm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * sin(m*lamda) ;
        dF_lamda_Snm =  (GM / l)   * (ae/l)^n * m * cos(m*lamda) * Pnm_norm(n+1,m+1) ; 
        
        % Cartesian counterparts of the partials
        fxyz_Cnm = PDVrx' * [dF_r_Cnm; dF_theta_Cnm; dF_lamda_Cnm];
        fxyz_Snm = PDVrx' * [dF_r_Snm; dF_theta_Snm; dF_lamda_Snm];
        
        i_coef = i_coef + 1;
        partials_c(:,i_coef) = fxyz_Cnm;
        partials_s(:,i_coef) = fxyz_Snm; 
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Differentiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 < 0
n
m
C20 = -4.841694947139e-04
Cnm(n+1,m+1) = C20;
S20 = 0;
Snm(n+1,m+1) = S20;
        dV_r = (GM/l^2) * ( -(n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) * (Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda)) );

delta_C = Cnm(n+1,m+1) * 10^-7
delta_S = 10^-7
        dV_r_dC = (GM/l^2) * ( -(n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) * ( (Cnm(n+1,m+1) + delta_C) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda)) );
        dV_r_dS = (GM/l^2) * ( -(n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) * ( (Cnm(n+1,m+1)) * cos(m*lamda) + ( Snm(n+1,m+1) + delta_S) * sin(m*lamda)) );

delta_dF_r_Cnm = (dV_r_dC - dV_r) / delta_C       
delta_dF_r_Snm = (dV_r_dS - dV_r) / delta_S  

test = stop 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

% Overall Partials matrix w.r.t. gravitational parameters
partials_param = [partials_c partials_s];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


% to_pod = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C coefficients
% i_coef = 0;
for i_C = 1 : Nparam_C
    n = C_degree_order(i_C, 1);
    m = C_degree_order(i_C, 2);
    
    dF_r_Cnm     = -(GM / l^2) * (n+1) * (ae/l)^n * Pnm_norm(n+1,m+1) * cos(m*lamda);
    dF_theta_Cnm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * cos(m*lamda) ;
    dF_lamda_Cnm =  (GM / l) * m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (- sin(m*lamda));    
    % Cartesian counterparts of the partials
    fxyz_Cnm = PDVrx' * [dF_r_Cnm; dF_theta_Cnm; dF_lamda_Cnm];
    
    % i_coef = i_coef + 1;
    % partials_param(:,i_coef) = fxyz_Cnm;
    
    partials_c(:,i_C) = fxyz_Cnm;
end

% S coefficients
for i_C = 1 : Nparam_S
    n = S_degree_order(i_C, 1);
    m = S_degree_order(i_C, 2);

    dF_r_Snm     = -(GM / l^2) * (n+1) * (ae/l)^n * Pnm_norm(n+1,m+1) * sin(m*lamda);
    dF_theta_Snm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * sin(m*lamda) ;
    dF_lamda_Snm =  (GM / l)   * (ae/l)^n * m * cos(m*lamda) * Pnm_norm(n+1,m+1) ;    
    % Cartesian counterparts of the partials
    fxyz_Snm = PDVrx' * [dF_r_Snm; dF_theta_Snm; dF_lamda_Snm];
    
    % i_coef = i_coef + 1;
    % partials_param(:,i_coef) = fxyz_Snm;
    partials_s(:,i_C) = fxyz_Snm;
end       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Overall matrix for C and S coefficients
partials_param = [partials_c partials_s]; 

% fprintf('%s %.3f \n', 'Time (sec): GRAV-partials 0:',toc(to_pod));


if 1 < 0
to_pod = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C coefficients
% i_coef = 0;
parfor i_C = 1 : Nparam_C
    n = C_degree_order(i_C, 1);
    m = C_degree_order(i_C, 2);
    
    dF_r_Cnm     = -(GM / l^2) * (n+1) * (ae/l)^n * Pnm_norm(n+1,m+1) * cos(m*lamda);
    dF_theta_Cnm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * cos(m*lamda) ;
    dF_lamda_Cnm =  (GM / l) * m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (- sin(m*lamda));    
    % Cartesian counterparts of the partials
    fxyz_Cnm = PDVrx' * [dF_r_Cnm; dF_theta_Cnm; dF_lamda_Cnm];
    
    % i_coef = i_coef + 1;
    % partials_param(:,i_coef) = fxyz_Cnm;
    
    partials_c(:,i_C) = fxyz_Cnm;
end

% S coefficients
parfor i_C = 1 : Nparam_S
    n = S_degree_order(i_C, 1);
    m = S_degree_order(i_C, 2);

    dF_r_Snm     = -(GM / l^2) * (n+1) * (ae/l)^n * Pnm_norm(n+1,m+1) * sin(m*lamda);
    dF_theta_Snm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * sin(m*lamda) ;
    dF_lamda_Snm =  (GM / l)   * (ae/l)^n * m * cos(m*lamda) * Pnm_norm(n+1,m+1) ;    
    % Cartesian counterparts of the partials
    fxyz_Snm = PDVrx' * [dF_r_Snm; dF_theta_Snm; dF_lamda_Snm];
    
    % i_coef = i_coef + 1;
    % partials_param(:,i_coef) = fxyz_Snm;
    partials_s(:,i_C) = fxyz_Snm;
end       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Overall matrix for C and S coefficients
partials_param_1 = [partials_c partials_s]; 

fprintf('%s %.3f \n', 'Time (sec): GRAV-partials 1:',toc(to_pod));


to_pod = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C and S coefficients
parfor i_C = 1 : Nparam_C
    n = C_degree_order(i_C, 1);
    m = C_degree_order(i_C, 2);
    
% C coefficients
    dF_r_Cnm     = -(GM / l^2) * (n+1) * (ae/l)^n * Pnm_norm(n+1,m+1) * cos(m*lamda);
    dF_theta_Cnm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * cos(m*lamda) ;
    dF_lamda_Cnm =  (GM / l) * m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (- sin(m*lamda));    
    % Cartesian counterparts of the partials
    fxyz_Cnm = PDVrx' * [dF_r_Cnm; dF_theta_Cnm; dF_lamda_Cnm];

% S coefficients
    dF_r_Snm     = -(GM / l^2) * (n+1) * (ae/l)^n * Pnm_norm(n+1,m+1) * sin(m*lamda);
    dF_theta_Snm =  (GM / l)   * (ae/l)^n * dPnm_norm(n+1,m+1) * sin(m*lamda) ;
    dF_lamda_Snm =  (GM / l)   * (ae/l)^n * m * cos(m*lamda) * Pnm_norm(n+1,m+1) ;    
    % Cartesian counterparts of the partials
    fxyz_Snm = PDVrx' * [dF_r_Snm; dF_theta_Snm; dF_lamda_Snm];
    
    partials_c(:,i_C) = fxyz_Cnm;
    partials_s(:,i_C) = fxyz_Snm;      
end   
% Overall matrix for C and S coefficients
partials_param_2 = [partials_c partials_s]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s %.3f \n', 'Time (sec): GRAV-partials 1:',toc(to_pod));

delta_partials = [partials_param_2 - partials_param_1]';
save delta_partials.out delta_partials -ASCII -double

partials_param = [partials_c partials_s]; 

end


%error
