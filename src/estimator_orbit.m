function [Xmatrix, Xmatrix_alt, bmatrix, Amatrix_out, Cx, Cv, NEQn, NEQu] = estimator_orbit (orbref,veqZarray,veqParray,orbobs,sigma_obs,obstype)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters Estimation
% Least-squares Estimator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - orbref:         Dynamic Orbit (observation function) 
%                   orbref = [t r' v'] 
%                   t: MJD incdluding fraction of the day 
%                   r,v: Position and Velcotiy vector cartesian coordinates 
% - veqZarray:      Variational Equations solution for the initial state
%                   vector (State transition matrix) per epoch                   
% - veqParray:      Variational Equations solution for additional parameters
%                   (force related parameters) per epoch 
% - orbobs:         Pseudo-Observations based on Kinematic Orbit positios (XYZ) 
%                   orbobs = [t r' v'] 
%                   t: MJD incdluding fraction of the day 
%                   r,v: Position and Velcotiy vector cartesian coordinates 
%
% Output arguments:
% - Xmatrix:        Estimated Parameters corrections 
% - Amatrix:        Design matrix
% - bmatrix:        b matrix of least squares method  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas D. Papanikolaou                                    August 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 15/04/2021, Thomas Papanikolaou
%             Upgrade and rename of function mainf_DOD.m
% 17/08/2022, Thomas Loudis Papanikolaou
%             Upgrade: matrix decomposition approaches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least-squares estimator
% Normal Equations are formed and computed directly at the common epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz1 sz2] = size(orbobs);
[sz3 sz4] = size(orbref);
iobs0 = 1;
Nepochs = 0;
for iref = 1 : sz3
    % Reference Orbit
    tiref = orbref(iref,1);
    % Pseduo-Observations (Kinematic Orbit)
    for iobs = iobs0 : sz1
        tiobs = orbobs(iobs,1);
        % Common Epochs
        if abs(tiref - tiobs) < 10^-8
            Nepochs = Nepochs + 1;
            % ObsEpochs(Nepochs,1) = orbref(iref,1);
            % obstype  :  1/2  [r/v]
            obstype = 1;
            if obstype == 1
                Nobsset = 3;
            elseif obstype == 2
                Nobsset = 6;
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Wmatrix_ti = [orbobs(iobs,2:(Nobsset+1)) - orbref(iref,2:(Nobsset+1))]';
            Wmatrix((Nepochs-1)*Nobsset+1 : Nepochs*Nobsset,1) = Wmatrix_ti;
            timearray = [tiref; tiref; tiref];
            Wmatrix_time((Nepochs-1)*Nobsset+1 : Nepochs*Nobsset,:) = [timearray Wmatrix_ti];           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Design Matrix values obtained by VEQ arrays values
            % veqZarray
            Amatrix_ti = veqZarray( (iref-1)*6+1 : (iref-1)*6+Nobsset , 2:end);
            Amatrix((Nepochs-1)*Nobsset+1 : Nepochs*Nobsset,:) = Amatrix_ti;
            Amatrix_time((Nepochs-1)*Nobsset+1 : Nepochs*Nobsset,:) = [timearray Amatrix_ti];
            % veqParray
            if size(veqParray,1) > 1 
            AmatrixP_ti = veqParray( (iref-1)*6+1 : (iref-1)*6+Nobsset , 2:end);
            AmatrixP((Nepochs-1)*Nobsset+1 : Nepochs*Nobsset,:) = AmatrixP_ti;
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            iobs0 = iobs + 1';
            clear tiref tiobs            
            break
        end
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design Matrix: Final Amatrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(veqParray,1) > 1 
Amatrix_ZP = [Amatrix AmatrixP];
else
Amatrix_ZP = Amatrix;
end
[d1,d2] = size(Amatrix_ZP);
Nobs = d1;
Nparam = d2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bmatrix = Wmatrix;
Amatrix_out = Amatrix_ZP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weighted Least-Squares method solution
[Xmatrix, NEQn, NEQu, error_matrix, sigma0, Cx, Cv] = estimator_neq_sol(Amatrix_out, bmatrix, sigma_obs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Errors matrix
% error_matrix = bmatrix - Amatrix_out * Xmatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariance matrix of estimated parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xmatrix_alt = Xmatrix;
