function [orbc,err,veqZarray,veqParray,orbk,orbt, forces_accel, Gmatrix, Rmatrix] = intg_main(intg,Zo,arc,MSprm,RKprm,eopdat,dpint,intgVEQ, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration of Equation of Motion and Variational Equations based on
% Multistep methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments
% - intg:         Intregrators ID
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
% - MSparam:      Multistep method's parameters
%   MSparam(1,1): Multistep method's ID number
%   MSparam(2,1): Order 
%   MSparam(3,1): Stepsize 
%   MSparam(4,1): Start integrator method for the first epochs
% - eopdat:       Earth Orientation Parameters (EOP) data that are required
%                 for the orbit arc length
% - dpint:        Number of data points (days) that are required for the 
%                 EOP interpolation to the computation epoch
% - intgVEQ:      Simple or Combined numerical integration
%   intgVEQ = 0 >> Integration of Equation of Motion
%   intgVEQ = 1 >> Integration of Equation of Motion and Variational
%   Equations in a simultaneous (combined) intgration formula
%
% Output arguments:
% - orbc:         State vector in the Celestial Reference System GCRS
%   orbc = [ t r_GCRS' v_GCRS' er' ev']
%   t:            Epoch in seconds of TT (Terrestrial Time)
% - err:          Integrator's local truncation errors  err = [er' ev']
%   er:           local truncation error of position
%   ev:           local truncation error of velocity
% - veqZarray:    VEQ array of state transition matrix (6*Epochs x 6)
% - veqParray:    VEQ array of sensitivity matrix      (6*Epochs x np)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                           May 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 07/12/2022, Dr. Thomas Loudis Papanikolaou
%             Code upgrade 
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Integration of Equation of Motion and Variational Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if intgVEQ == 0
    % Equation of motion
    % (Single numerical integration)
    if intg == 1 || intg == 2 || intg == 3
        [orbc,err] = integr_rk(Zo,arc,RKprm,eopdat,dpint, orbit_model_struct);
        forces_accel = 0;
    else
        [orbc, forces_accel] = integr_ms(Zo,arc,MSprm,RKprm,eopdat,dpint,orbit_model_struct);
        err = 0;
    end
    veqZarray = 0;
    veqParray = 0;
elseif intgVEQ == 1
    % Equation of motion and Variational Equations
    % (Combined numerical integration) 
    if intg == 1 || intg == 2 || intg == 3
        [orbc,err,veqZarray,veqParray] = integr_rk_veq(Zo,arc,RKprm,eopdat,dpint, orbit_model_struct);
    else
        [orbc,err,veqZarray,veqParray] = integr_ms_veq(Zo,arc,MSprm,RKprm,eopdat,dpint, orbit_model_struct);
    end
    forces_accel = 0;
end

if MSprm(1,1) == 6 || MSprm(1,1) == 7
    %fprintf('%s \n', 'Stoermer-Cowell methods')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keplerian Elements computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force model structure matrix
GM_glob = orbit_model_struct.GM_Earth;

[sz1 sz2] = size(orbc);
% Preallocation 
orbk = zeros(sz1, sz2);
% Keplerian Elements array
for ik = 1 : sz1
    [a, e, incl, Omega, omega, f, Mn, E, u] = kepler(orbc(ik,2:4)',orbc(ik,5:7)',GM_glob);
    orbk(ik,:) = [orbc(ik,1) a e incl Omega omega Mn];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit transformation from GCRS to ITRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[orbt] = orbc2t(orbc,eopdat,dpint,orbit_model_struct);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Gradient matrix and Earth Orientation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gmatrix = zeros(3,3);
Rmatrix = zeros(3,3);

