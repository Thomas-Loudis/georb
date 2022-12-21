function [orbc,err] = integr_rk(zo,arc,RKparam,eop,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration of Equation of Motion based on Runge-Kutta-Nystrom methods
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments
% - zo:           Initial epoch and state vector in Celestial Reference
%                 System GCRS
%   zo = [to ro' vo']
%   to:           MJD of initial epoch in days
%   ro:           initial position in GCRS (m)
%   vo:           initial velocity in GCRS (m/sec)
% - arc:          Orbit arc lenth in seconds
% - RKparam:      Runge-Kutta & Runge-Kutta-Nystrom methods parameters
%   RKparam(1,1): RK method's ID number
%   RKparam(2,1): Stepsize (h)
%   RKparam(3,1): RKN7(6)-8 parameter for stepsize control (lamda)
%   RKparam(4,1): RKN6(4)-4FD parameter for interpolation points
%                 RKparam(4,1) = sigma, 0<sigma<h
% - eop:          Earth Orientation Parameters (EOP) data that are required
%                 for the orbit arc length
% - dpint:        Number of data points (days) that are required for the 
%                 EOP interpolation to the computation epoch
%
% Output arguments:
% - orbc:         Array of MJD, State vector, Integrator errors in the GCRS
%   orbc = [t r_GCRS' v_GCRS' er' ev']
%   t:            MJD of epochs in days
% - err:          Integrator's local truncation errors  err = [er' ev']
%   er:           local truncation error of position
%   ev:           local truncation error of velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   June  2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
%   06/05/2010    Minor upgrades. Function has been renamed.
%   25/04/2011    Compatible with the "new parameterization", the new
%                 functions for IERS-EOP and especially the time input
%                 argument based on MJD.          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RK integration method's parameters
h = RKparam(2,1);
sigma = RKparam(4,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time conversion to seconds (MJD from days to seconds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial epoch
to = zo(1,1) * (24*3600);
tmax = to + arc;
Nepochs = fix(arc/h)+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation array
orbc = zeros(Nepochs,13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial position and velocity vectors in the GCRS
ro = [zo(1,2); zo(1,3); zo(1,4)];
vo = [zo(1,5); zo(1,6); zo(1,7)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial State Vector in the GCRS
zGCRS = [to ro' vo'];
orbc(1,:) = [zGCRS 0 0 0 0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration method is defined by the element "RKparam(1,1)"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RKN7(6)-8 method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RKparam(1,1) == 2
i = 1;
TT = to;
for i = 2 : Nepochs
    % computation of state vector at next epoch in the GCRS
    [z_q,e_z] = rkn768(zGCRS,RKparam,eop,dpint);
    % State vector at next epoch TT (to+h) in the GCRS
    TT = TT + h;    
    rGCRS = z_q(1,1:3)';
    vGCRS = z_q(1,4:6)';
    zGCRS = [TT rGCRS' vGCRS'];
    orbc(i,:) = [zGCRS e_z];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RKN6(4)-4FD  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif RKparam(1,1) == 3
i = 1;
for t = to : h : tmax-1
    i = i + 1;
    % computation of state vector at next epoch in the GCRS
    [z_q,e_z,z_int,ez_int] = rkn646fd(zGCRS,RKparam,eop,dpint);
    % State vector at next epoch (to+h) in the GCRS
    TT = t + h;
    rGCRS = z_q(1,1:3)';
    vGCRS = z_q(1,4:6)';
    zGCRS = [TT rGCRS' vGCRS'];    
    orbc(i,:) = [zGCRS e_z];      
    % RKN method interpolant
    %  State vector at interpolation epoch (to+sigma*h) in the GCRS
    TT_int = t + sigma * h;
    rGCRS_int = z_int(1,1:3)';
    vGCRS_int = z_int(1,4:6)';
    orbc_int(i,:) = [TT_int rGCRS_int' vGCRS_int' ez_int];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epochs MJD conversion from seconds to days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orbc(:,1) = orbc(:,1) / (24 * 3600);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err = orbc(:,8:13);
orbc = orbc(:,1:7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
