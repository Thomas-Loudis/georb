function [z_q,e_z] = rkn768(zo,RKparam,eop,dpint,orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta-Nystrom methods : 
% Implementation of RKN7(6)-8 method by Dormand & Prince(1987) for the
% solution of the Equation of Motion 
%
% RKN methods are especially designed for the direct integration of
% second-order differential equation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
%  RKN7(6)-8 method:
%  Accelerometry does not depend on the velocity of the celestial body.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References
%  Dormand J.R., Prince P.J. (1987). Runge-Kutta-Nystrom triples. Comput.
%  Math. Applic., Vol.13, No.12, pp. 937-949.
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
% Output arguments
% - z_q:   State vector in GCRS at next epoch t=to+h
%          z_q = [t r_q' v_q']
%   r_q:   position at epoch t=to+h in GCRS
%   v_q:   velocity at epoch t=to+h in GCRS
% - e_z:   Local truncation error
%          e_z = [e_r' e_v']
%   er:    local truncation error of position
%   ev:    local truncation error of velocity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
% Function 'accel.m' is called for the computation of the acceleration in
% the Inertial System (GCRS).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                     September 2007
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
for i = 0 : s
    if i == 0
%      Function evaluations (ki) for State vector (z = [r v]')       
        t(i+1) = to;        
        r(:,i+1) = ro;
        % random value for velocity v = vo:
        % Force model is indepedent from velocity vector
        zmjd = [t(i+1)/(24*3600) r(:,i+1)' vo'];
        [fx,fy,fz] = accel(zmjd,eop,dpint,orbit_model_struct);
        k(:,i+1) = [fx; fy; fz];
    else        
%      Function evaluations (ki) for State vector        
        for j = 0 : i-1
            sum_ak = sum_ak + a(i+1,j+1) * k(:,j+1);
        end        
        t(i+1) = to + c(i+1,1) * h;
        r(:,i+1) = ro + c(i+1,1) * h * vo + h^2 * sum_ak ;
        sum_ak = 0;
        % random value for velocity v = vo:
        % Force model is indepedent from velocity vector
        zmjd = [t(i+1)/(24*3600) r(:,i+1)' vo'];
        [fx,fy,fz] = accel(zmjd,eop,dpint, orbit_model_struct);
        k(:,i+1) = [fx; fy; fz];
    end
end
% Increment function for State vector
sum_bk = [0;0;0] ;
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