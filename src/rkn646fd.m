function [z_q,e_z,z_int,ez_int] = rkn646fd(zo,RKparam,eop,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runge-Kutta-Nystrom methods
% 
% Implementation of RKN6(4)-6FD method (Continuous method - Interpolant) by
% Dormand & Prince(1978) for the solution of the Equation of Motion 
%
% RKN methods are especially designed for the direct integration of
% second-order differential equation
%
% RKN6(4)-4FD method:
%  Accelerometry does not depend on the velocity of the celestial body
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References
%  Dormand J.R., Prince P.J. (1978). New Runge-Kutta algorithms for
%  numerical simulation in dynamical astronomy. Celestial Mechanics, 18,
%  pp. 223-232.   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments
% - zo:  Initial epoch and state vector in Celestial Reference System GCRS
%        zo = [to ro' vo']  
%   to:  initial epoch in seconds of TT (Terrestrial Time)
%   ro:  initial position in GCRS (m)
%   vo:  initial velocity in GCRS (m/sec)
%
% Output arguments
% - z_q:     State vector in GCRS at next epoch t=to+h
%            z_q = [t r_q' v_q']
%   r_q:     Position at epoch t=to+h in GCRS
%   v_q:     Velocity at epoch t=to+h in GCRS
% - e_z:     Local truncation error
%            e_z = [e_r' e_v']
%   er:      local truncation error of position
%   ev:      local truncation error of velocity
% - z_int:   State vector in GCRS at interpolation points (t=sigma*h)
%   r_int:   position at interpolation epochs in GCRS
%   v_int:   velocity at interpolation epochs in GCRS
% - ez_int:  Local truncation error at interpolation points
%            ez_int = [er_int' ev_int'];
%   er_int:  local truncation error of position
%   ev_int:  local truncation error of velocity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remarks:
% Function m.files that are called during computations
% - 'Acceleration.m' is called for the computation of the acceleration
%   in the Inertial System (GCRS).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                     November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RK integration method's parameters
h = RKparam(2,1);
sigma = RKparam(4,1);
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
% s-stage function evaluations k(s), s=0..to...
s = 6;
R = sqrt(8581);
% ci coefficients, i = 1,....,s
c=[ 0
    (209-R)/900
    (209-R)/450
    (209+R)/450
    9/10
    1];
% aij coefficients, i=2,....,s, j=1,...,s-1
a=[                   NaN                                            NaN                                                   NaN                                              NaN                            NaN       
             (26131-209*R)/810000                                    NaN                                                   NaN                                              NaN                            NaN
             (26131-209*R)/607500                            (26131-209*R)/303750                                          NaN                                              NaN                            NaN
    (980403512254+7781688431*R)/11694469921875     -(1262884486208+15385481287*R)/11694469921875      (7166233891441+78694563299*R)/46777879687500                          NaN                            NaN
          -9*(329260+3181*R)/27040000                      27*(35129+3331*R)/13520000                    -27*(554358343+31040327*R)/464060480000            153*(8555257-67973*R)/2745920000               NaN
                   329/4212                                           0                                       (84119543+366727*R)/409622616                   (84119543-366727*R)/409622616             200/17901];
% position coefficients bi, order q (q>p),  i=1,....,s
bq_r = [329/4212 0 (84119543+366727*R)/409622616 (84119543-366727*R)/409622616 200/17901 0]';
% velocity coefficients bi, order q (q>p),  i=1,....,s
bq_v = [329/4212 0 (389225579+96856*R)/1024056540 (389225579-96856*R)/1024056540 2000/17901 1/20]';
% position coefficients bi, order p,  i=1,....,s
bp_r = [(2701+23*R)/4563 -(9829+131*R)/9126 5*(1798+17*R)/9126 0 0 0 ]';
% velocity coefficients bi, order p,  i=1,....,s
bp_v = [115/2106 0 (84119543+366727*R)/256014135 (84119543-366727*R)/256014135 6950/17901 -1/10 ]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolant coefficients 
bint_r=[(900*sigma^4-3819*sigma^3+6386*sigma^2-5244*sigma+2106)/4212
        0
        sigma*(1800*(5860823-152228*R)*sigma^3-6*(4929647204-156109769*R)*sigma^2+(22190560391-1109665151*R)*sigma+18*(81356461+25954829*R))/22529243880 
        sigma*(1800*(5860823+152228*R)*sigma^3-6*(4929647204+156109769*R)*sigma^2+(22190560391+1109665151*R)*sigma+18*(81356461-25954829*R))/22529243880 
       ( -200*sigma*(225*sigma^3-651*sigma^2+620*sigma-195) )/17901
       (sigma*(sigma-1)*(300*sigma^2-523*sigma+234))/220];

bint_v=[(5400*sigma^4-19095*sigma^3+25544*sigma^2-15732*sigma+4212)/4212
         0
         sigma*(5400*(5860823-152228*R)*sigma^3-15*(4929647204-156109769*R)*sigma^2+2*(22190560391-1109665151*R)*sigma+27*(81356461+25954829*R))/11264621940
         sigma*(5400*(5860823+152228*R)*sigma^3-15*(4929647204+156109769*R)*sigma^2+2*(22190560391+1109665151*R)*sigma+27*(81356461-25954829*R))/11264621940
         (-1000*sigma*(270*sigma^3-651*sigma^2+496*sigma-117))/17901
         (sigma*(1800*sigma^3-4115*sigma^2+3028*sigma-702))/220]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations for the description of the Runge-Kutta-Nystrom method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function evaluations K(i), i=1,...,s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K(i) is expressed here as matrix K(:,i) = [Kx(i); Ky(i); Kz(i)]
sum_ak = [0;0;0];
for i = 1 : s
    if i == 1
%      Function evaluations (ki) for State vector (z = [r v]')       
        t(i) = to;        
        r(:,i) = ro;
        % random value for velocity v = vo:
        % Force model is indepedent from velocity vector 
        %[fx,fy,fz] = Acceleration(t(i),r(:,i),vo);
        zmjd = [t(i)/(24*3600) r(:,i)' vo'];
        [fx,fy,fz] = accel(zmjd,eop,dpint);
        k(:,i) = [fx; fy; fz];
        %k(:,i)= [0; 0; 0;];       
    else
%      Function evaluations (ki) for State vector                
        for j=1:i-1           
            sum_ak = sum_ak + a(i,j) * k(:,j);
        end
        t(i) = to + c(i,1) * h;
        r(:,i) = ro + c(i,1) * h * vo + h^2 * sum_ak;        
        sum_ak = 0;        
        % random value for velocity v = vo:
        % Force model is indepedent from velocity vector 
        %[fx,fy,fz] = Acceleration(t(i),r(:,i),vo);
        zmjd = [t(i)/(24*3600) r(:,i)' vo'];
        [fx,fy,fz] = accel(zmjd,eop,dpint);        
        k(:,i) = [fx; fy; fz];
    end
end
% Increment function for State vector
sum_bq_r_k = [0;0;0];
sum_bq_v_k = [0;0;0];
sum_bp_r_k = [0;0;0];
sum_bp_v_k = [0;0;0];
for i = 1 : s
    sum_bq_r_k = sum_bq_r_k + bq_r(i,1) * k(:,i);
    sum_bq_v_k = sum_bq_v_k + bq_v(i,1) * k(:,i);
    sum_bp_r_k = sum_bp_r_k + bp_r(i,1) * k(:,i);
    sum_bp_v_k = sum_bp_v_k + bp_v(i,1) * k(:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution
%  Embedded Runge-Kutta method
%  Position and velocity approximations of order q
r_q = ro + h * vo + (h)^2 * sum_bq_r_k;
v_q = vo + h * sum_bq_v_k;
z_q = [r_q' v_q'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local truncation error
%  Position approximations of lower order p
r_p = ro + h * vo + (h)^2*sum_bp_r_k;
v_p = vo + h * sum_bp_v_k;
%  Local truncation error - Stepsize control
er = abs(r_q - r_p);
ev = abs(v_q - v_p);
e_z = [er' ev'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dense formula
% dense output points are defined by step-size sigma*h (0<sigma<1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remarks:
%  RKN6(4)-4FD method's interpolant function evaluations are equal with
%  the evaluations that have been computed during the integration step.
%  Thus, no extra function evaluations are required for the interpolant.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_bint_r_k = [0; 0; 0];
sum_bint_v_k = [0; 0; 0];
for i = 1 : s
    sum_bint_r_k = sum_bint_r_k + bint_r(i,1)*k(:,i);
    sum_bint_v_k = sum_bint_v_k + bint_v(i,1)*k(:,i);
end
% Position and velocity approximations of order q at interpolation points
r_int = ro + sigma *h*vo + (sigma * h)^2 * sum_bint_r_k ;
v_int = vo + sigma *h*sum_bint_v_k ;
z_int = [r_int' v_int'];
% Local truncation errors at interpolation points
er_int = [0; 0; 0];
ev_int = [0; 0; 0];
ez_int = [er_int' ev_int'];