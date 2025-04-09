function [z_q,e_z,veqZ,veqP] = rkn768_veq(zo,RKparam,eop,dpint,veqZo,veqPo, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined Integration of Equation of Motion and Variational Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta-Nystrom methods:
% Implementation of RKN7(6)-8 method by Dormand & Prince(1987) for the
% solution of the Equation of Motion and Variational Equations
%
% RKN methods are especially designed for the direct integration of
% second-order differential equation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
% In RKN methods, acceleration does not depend on the velocity of the
% celestial body.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References
%  Dormand J.R., Prince P.J. (1987). Runge-Kutta-Nystrom triples. Comput.
%  Math. Applic., Vol.13, No.12, pp. 937-949
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments
% - zo:  Initial epoch and state vector in Celestial Reference System GCRS
%        zo = [to ro' vo']  
%   to:  initial epoch in seconds of TT (Terrestrial Time)
%   ro:  initial position in GCRS (m)
%   vo:  initial velocity in GCRS (m/sec)
% - RKparam:       Runge-Kutta & Runge-Kutta-Nystrom methods parameters
%   RKparam(1,1):  RK method's ID number
%   RKparam(2,1):  Stepsize (h)
%   RKparam(3,1):  RKN7(6)-8 parameter for stepsize control (lamda)
%   RKparam(4,1):  RKN6(4)-4FD parameter for interpolation points
%                  RKparam(4,1) = sigma, 0<sigma<h
% - eop:    Earth Orientation Parameters (EOP) data that are required for
%           the orbit arc length
% - dpint:  Number of data points (days) that are required for the EOP
%           interpolation to the computation epoch
% - veqZo:  State transition matrix at initial epoch (6x6 matrix)
% - veqPo:  Sensitivity matrix at initial epoch (6xnp matrix)
%
% Output arguments
% - z_q:   State vector in GCRS at next epoch t=to+h
%          z_q = [t r_q' v_q']
%   r_q:   position at epoch t=to+h in GCRS
%   v_q:   velocity at epoch t=to+h in GCRS
% - e_z:   Local truncation error
%          e_z = [e_r' e_v']
%   er:    local truncation error of position
%   ev:    local truncation error of velocity
% - veqZ:  State transition matrix at computation epoch (6x6 matrix)
% - veqP:  Sensitivity matrix at computation epoch (6xnp matrix)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State transition matrix: STM matrix 6x6
% 
% STM = [ STMxxo  STMxyo  STMxzo  STMxVxo  STMxVyo  STMxVzo
%         STMyxo  STMyyo  STMyzo  STMyVxo  STMyVyo  STMyVzo
%         STMzxo  STMzyo  STMzzo  STMzVxo  STMzVyo  STMzVzo  
%         STMVxxo STMVxyo STMVxzo STMVxVxo STMVxVyo STMVxVzo  
%         STMVyxo STMVyyo STMVyzo STMVyVxo STMVyVyo STMVyVzo  
%         STMVzxo STMVzyo STMVzzo STMVzVxo STMVzVyo STMVzVzo ]
%
% STM_r = [ STMxxo  STMxyo  STMxzo  STMxVxo  STMxVyo  STMxVzo
%           STMyxo  STMyyo  STMyzo  STMyVxo  STMyVyo  STMyVzo
%           STMzxo  STMzyo  STMzzo  STMzVxo  STMzVyo  STMzVzo ]
%
% STM_v = [ STMVxxo STMVxyo STMVxzo STMVxVxo STMVxVyo STMVxVzo  
%           STMVyxo STMVyyo STMVyzo STMVyVxo STMVyVyo STMVyVzo  
%           STMVzxo STMVzyo STMVzzo STMVzVxo STMVzVyo STMVzVzo ]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                     September 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 07/05/2012  Upgrade for (re)development of orbit adjustment algorithm
% 04/06/2012  Replacement of the functions accel.m and pdv_accl.m by the
%             combined function veq_accl.m
% 22/08/2013  Upgrade for adding empirical modelling parameters (e.g. 1CPR) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RK integration method's parameters
h = RKparam(2,1);
lamda_h = RKparam(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial epoch
to = zo(1,1);
% Initial position and velocity vectors
ro = zo(1,2:4)';
vo = zo(1,5:7)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State Transition Matrix Initial values: veqZo = I (6x6)
veq_ro = zeros(3,6);
veq_vo = zeros(3,6);
ir = 1;
iv = 1;
for i1 = 1 : 3
    veq_ro(ir,:) = veqZo(i1,:);
    ir = ir + 1;
end
clear i1 ir
for i1 = 4 : 6
    veq_vo(iv,:) = veqZo(i1,:);
    iv = iv + 1;
end
clear i1 iv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity Matrix Initial values: veqPo (6xNp)
[sz6, Np_veqP] = size(veqPo);
veqP_ro = zeros(3,Np_veqP);
veqP_vo = zeros(3,Np_veqP);
ir = 1;
for i1 = 1 : 3
    veqP_ro(ir,:) = veqPo(i1,:);
    ir = ir + 1;
end
clear i1 ir
iv = 1;
for i1 = 4 : 6
    veqP_vo(iv,:) = veqPo(i1,:);
    iv = iv + 1;
end
clear i1 iv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nparam = Np_veqP;
% Matrices prealocation
gZ_3d = zeros(3,6,9);
gP_3d = zeros(3,Nparam,9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients ci, aij, bi 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s-stage Function evaluations k(s), s=0..to...
s = 8;
% order p=7

% ci coefficients, i=0,....,s
c=[ 0
    1/10
    1/5
    3/8
    1/2
    (7-sqrt(21))/14
    (7+sqrt(21))/14
    1
    1];
% aij coefficients, i=1,....,s, j=0,...,s-1
a=[  NaN                         NaN                            NaN                             NaN                                 NaN                          NaN                          NaN                  NaN
    1/200                        NaN                            NaN                             NaN                                 NaN                          NaN                          NaN                  NaN      
    1/150                        1/75                           NaN                             NaN                                 NaN                          NaN                          NaN                  NaN 
    171/8192                     45/4096                        315/8192                        NaN                                 NaN                          NaN                          NaN                  NaN      
    5/288                        25/528                         25/672                          16/693                              NaN                          NaN                          NaN                  NaN 
    (1003-205*sqrt(21))/12348   -25*(751-173*sqrt(21))/90552    25*(624-137*sqrt(21))/43218    -128*(361-79*sqrt(21))/237699        (3411-745*sqrt(21))/24696    NaN                          NaN                  NaN 
    (793+187*sqrt(21))/12348    -25*(331+113*sqrt(21))/90552    25*(1044+247*sqrt(21))/43218   -128*(14885+3779*sqrt(21))/9745659   (3327+797*sqrt(21))/24696   -(581+127*sqrt(21))/1722      NaN                  NaN  
    -(157-3*sqrt(21))/378        25*(143-10*sqrt(21))/2772     -25*(876+55*sqrt(21))/3969       1280*(913+18*sqrt(21))/596673      -(1353+26*sqrt(21))/2268      7*(1777+377*sqrt(21))/4428   7*(5-sqrt(21))/36    NaN
    1/20                         0                              0                               0                                   8/45                         7*(7+sqrt(21))/360           7*(7-sqrt(21))/360   0   ];
% position coefficients bi, order p+1,  i=0,....,s
b1 = [ 1/20 0 0 0 8/45 7*(7+sqrt(21))/360 7*(7-sqrt(21))/360 0 0 ]';
% position coefficients bi, order p,  i=0,....,s
b = [ 1/20 0 0 0 8/45 7*(7+sqrt(21))/360 7*(7-sqrt(21))/360 -lamda_h lamda_h ]';
% velocity coefficients bi, order p+1(p),  i=0,....,s
b2 = [ 1/20 0 0 0 16/45 49/180 49/180 1/20 0 ]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations for the description of the Runge-Kutta-Nystrom method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function evaluations K(i), i = 0,...,s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K(i) is expressed here as matrix K(:,i) = [Kx(i); Ky(i); Kz(i)]
k = zeros(3,s+1);
t = zeros(s+1);
r = zeros(3,s+1);
sum_ak = [0;0;0];
sum_ag = zeros(3,6);
sum_agP = zeros(3,Np_veqP); 
for i = 0 : s
    if i == 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Function evaluations (ki) for State vector (z = [r v]')       
        t(i+1) = to;        
        r(:,i+1) = ro;
        % random value for velocity v = vo:
        % Force model is indepedent from velocity vector
        zmjd = [t(i+1)/(24*3600) r(:,i+1)' vo'];
        % Acceleration & Partial Derivatives
        [fx,fy,fz,pdv_acc,pdv_acc_param] = veq_accl(zmjd,eop,dpint, orbit_model_struct);
        k(:,i+1) = [fx; fy; fz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Function evaluations (gi) for State Transition matrix
        %[pdv_acc] = pdv_accl(zmjd,eop,dpint);
        % 2nd derivative of veq matrices (fundamental VEQ formula)
        drv2_veq_r = pdv_acc * veq_ro;  % 3x6 matrix
        g = drv2_veq_r;
                
        % Function evaluation as a 3x6xi matrix
		% k_Z(:,:,i+1) = kZ_i(:,:)  										 
        gZ_3d(:,:,i+1) = g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Function evaluations (pi) for Sensitivity matrix
        %[pdv_acc] = pdv_accl(zmjd,eop,dpint);
        pdv_acc_P = pdv_acc_param;

        % 2nd derivative of veq matrices (fundamental VEQ formula)
        drv2_veqP_r = pdv_acc * veqP_ro + pdv_acc_P;  % 3xNp matrix
        gp = drv2_veqP_r;
               
        % Function evaluation as a 3xPxi matrix
		% k_P(:,:,i+1) = kP_i(:,:) 									
        gP_3d(:,:,i+1) = gp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Function evaluations (ki) for State vector        
        for j = 0 : i-1
            sum_ak = sum_ak + a(i+1,j+1) * k(:,j+1);
        end
        clear j
        t(i+1) = to + c(i+1,1) * h;
        r(:,i+1) = ro + c(i+1,1) * h * vo + h^2 * sum_ak ;
        
        % Compute velocity at intermediate epoch ti
        vi = vo + (h / c(i+1,1)) * sum_ak ;

        sum_ak = [0;0;0];

        % Force model is indepedent from velocity vector
        zmjd = [t(i+1)/(24*3600) r(:,i+1)' vi'];
        
        % Acceleration & Partial Derivatives
        [fx,fy,fz,pdv_acc,pdv_acc_param] = veq_accl(zmjd,eop,dpint, orbit_model_struct);
        k(:,i+1) = [fx; fy; fz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Function evaluations (gi) for State Transition matrix
        for j=0:i-1
            G_gZ = gZ_3d(:,:,j+1);     
            sum_ag = sum_ag + a(i+1,j+1) * G_gZ;
        end
        j=0;
    
        veq_r = veq_ro + c(i+1,1) * h * veq_vo + h^2 * sum_ag ;
        sum_ag = zeros(3,6);

        % 2nd derivative of veq matrices (fundamental VEQ formula)
        %[pdv_acc] = pdv_accl(zmjd,eop,dpint);
        drv2_veq_r = pdv_acc * veq_r;  % 3x6 matrix       
        g = drv2_veq_r;       
                
        % Function evaluation as a 3x6xi matrix
		% k_Z(:,:,i+1) = kZ_i(:,:)  										 
        gZ_3d(:,:,i+1) = g;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Function evaluations (pi) for Sensitivity matrix        
        for j=0:i-1
            G_gP = gP_3d(:,:,j+1);     
            sum_agP = sum_agP + a(i+1,j+1) * G_gP;
        end
        j=0;
    
        veqP_r = veqP_ro + c(i+1,1) * h * veqP_vo + h^2 * sum_agP ;
        sum_agP = zeros(3,Np_veqP);

        % 2nd derivative of veq matrices (fundamental VEQ formula)
        %[pdv_acc] = pdv_accl(zmjd,eop,dpint);        
        pdv_acc_P = pdv_acc_param; 
        
        drv2_veqP_r = pdv_acc * veqP_r + pdv_acc_P;  % 3xNp matrix       
        gp = drv2_veqP_r;    
        
        % Function evaluation as a 3xPxi matrix
		% k_P(:,:,i+1) = kP_i(:,:) 									
        gP_3d(:,:,i+1) = gp;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Increment function for State vector
sum_bk =  [0;0;0] ;
sum_b1k = [0;0;0] ;
sum_b2k = [0;0;0] ;
for i = 0 : s        
    sum_bk = sum_bk + b(i+1,1) * k(:,i+1);
    sum_b1k = sum_b1k + b1(i+1,1) * k(:,i+1);
    sum_b2k = sum_b2k + b2(i+1,1) * k(:,i+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution
%  Embedded Runge-Kutta method
%  Position and velocity approximations of order q(p+1)
r_q = ro + h * vo + (h)^2 * sum_b1k;
v_q = vo + h * sum_b2k;
z_q = [r_q' v_q'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local truncation error
%  Position approximation of lower order p
r_p = ro + h * vo + (h)^2 * sum_bk;
%  Local truncation error - Stepsize control
er = abs(r_q - r_p);
ev = [0; 0; 0];
e_z = [er' ev'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State Transition matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Increment function for State Transition matrix
sum_bq_r_g = [ zeros(3) zeros(3) ];
sum_bq_v_g = [ zeros(3) zeros(3) ];
sum_bp_r_g = [ zeros(3) zeros(3) ];
sum_bp_v_g = [ zeros(3) zeros(3) ];
    
for i = 0 : s       
    G_gZ = gZ_3d(:,:,i+1);        
    sum_bq_r_g = sum_bq_r_g + b1(i+1,1) * G_gZ;
    sum_bq_v_g = sum_bq_v_g + b2(i+1,1) * G_gZ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% approximations of State transition matrix of order q
veq_r = veq_ro + h * veq_vo + (h)^2 * sum_bq_r_g;
veq_v = veq_vo + h * sum_bq_v_g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
veqZ = [veq_r; veq_v];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Increment function for Sensitivity matrix
sum_bq_r_gP = [ zeros(3,Nparam) ];
sum_bq_v_gP = [ zeros(3,Nparam) ];
sum_bp_r_gP = [ zeros(3,Nparam) ];
sum_bp_v_gP = [ zeros(3,Nparam) ];
for i = 0 : s
    G_gP = gP_3d(:,:,i+1);
    sum_bq_r_gP = sum_bq_r_gP + b1(i+1,1) * G_gP;
    sum_bq_v_gP = sum_bq_v_gP + b2(i+1,1) * G_gP;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% approximations of State transition matrix of order q
veqP_r = veqP_ro + h * veqP_vo + (h)^2 * sum_bq_r_gP;
veqP_v = veqP_vo + h * sum_bq_v_gP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
veqP = [veqP_r; veqP_v];
