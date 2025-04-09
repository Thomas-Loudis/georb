function [orbc,err,veqZarray,veqParray] = integr_rk_veq(zo,arc,RKparam,eop,dpint, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration of Equation of Motion and Variational Equations based on
% Runge-Kutta-Nystrom methods
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
% - veqZarray:    VEQ array of state transition matrix (6*Epochs x 6)
% - veqParray:    VEQ array of sensitivity matrix      (6*Epochs x np)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                    September 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
%  07/05/2012  Thomas Papanikolaou
%              Upgrade for (re)development of orbit estimator algorithm
%  22/08/2013  Dr. Thomas Papanikolaou
%              Upgrade for adding empirical modelling parameters (e.g. CPR) 
%  22/08/2013  Dr. Thomas Papanikolaou
%              Upgrade for adding empirical modelling parameters (e.g. CPR) 
%  26/05/2020  Dr. Thomas Papanikolaou
%              Call variable Nparam_GLOB of the number of unkown parametrs
%               to be estimated in addition to the initial state vector             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nparam = orbit_model_struct.forces_param_estim_no;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation array
Nepochs = fix(arc/h)+1;
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
% VEQ initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation VEQ arrays
veqZarray = zeros(6*Nepochs, 1 + 6);
veqParray = zeros(6*Nepochs, 1 + Nparam);             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State Transition Matrix Initial values: veqZo = I (6x6)
veq_ro = [ 1 0 0 0 0 0
           0 1 0 0 0 0
           0 0 1 0 0 0 ];
veq_vo = [ 0 0 0 1 0 0
           0 0 0 0 1 0
           0 0 0 0 0 1 ];

veqZo = [veq_ro; veq_vo];

% Sensitivity Matrix initial values: veqPo = 0 (6xNp)
veqPo = zeros(6,Nparam);                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t6x1 = to * [1; 1; 1; 1; 1; 1];
veqZarray(1:6,:) = [t6x1 veqZo];
veqParray(1:6,:) = [t6x1 veqPo];
clear t6x1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RKN7(6)-8 method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RKparam(1,1) == 2
i = 1;
for t = to : h : tmax-1
    i = i + 1;        
    % computation of state vector at next epoch in the GCRS
    [z_q,e_z,veqZ,veqP] = rkn768_veq(zGCRS,RKparam,eop,dpint,veqZo,veqPo, orbit_model_struct);    
    % State vector at next epoch TT (to+h) in the GCRS
    TT = t + h;    
    rGCRS = z_q(1,1:3)';
    vGCRS = z_q(1,4:6)';
    zGCRS = [TT rGCRS' vGCRS'];
    orbc(i,:) = [zGCRS e_z];
    % VEQ arrays    
    t6x1 = TT * [1; 1; 1; 1; 1; 1];
    veqZarray( (i-1)*6+1 : (i-1)*6+6 , : ) = [t6x1 veqZ];
    veqParray( (i-1)*6+1 : (i-1)*6+6 , : ) = [t6x1 veqP];
    % Initial VEQ arrays at next epoch
    veqZo = veqZ;
    veqPo = veqP;
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
    [z_q,e_z,z_int,ez_int] = rkn646fd(zGCRS,RKparam,eop,dpint, orbit_model_struct);
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
veqZarray(:,1) = veqZarray(:,1) / (24 * 3600);
veqParray(:,1) = veqParray(:,1) / (24 * 3600);
