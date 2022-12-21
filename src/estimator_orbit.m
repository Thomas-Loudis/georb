function [Xmatrix,Xmatrix_alt, bmatrix, Amatrix, Cx, Cv] = estimator_orbit (orbref,veqZarray,veqParray,orbobs,COVobs,obstype)


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
            ObsEpochs(Nepochs,1) = orbref(iref,1);
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
            AmatrixP_ti = veqParray( (iref-1)*6+1 : (iref-1)*6+Nobsset , 2:end);
            AmatrixP((Nepochs-1)*Nobsset+1 : Nepochs*Nobsset,:) = AmatrixP_ti;
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
Amatrix_ZP = [Amatrix AmatrixP];
[d1,d2] = size(Amatrix_ZP);
Nobs = d1;
Nparam = d2;
Amatrix_all   = Amatrix_ZP;
Wmatrix_final = Wmatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bmatrix = Wmatrix_final;
Amatrix = Amatrix_all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least-Squares solution
[Xmatrix NEQn, NEQu] = estimator_neq_sol(Amatrix, bmatrix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Errors matrix
error_matrix = bmatrix - Amatrix * Xmatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariance matrix of estimated parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Observations
[dim1 dim2] = size(bmatrix);
n_obs = dim1;

% Number of Parameters 
[dim1 dim2] = size(Xmatrix);
m_param = dim1;

% Reference Variance
sigma = sqrt(error_matrix' * error_matrix / (n_obs - m_param) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariance matrix of parameters
%Cx = sigma^2 * inv(NEQn);

% Covariance matrix of errors
%Cv = sigma^2 * (inv(P_matrix) - Amatrix * inv(NEQn) * Amatrix');
%Cv = sigma^2 * ( - Amatrix * inv(NEQn) * Amatrix');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cx = 0;
Cv = 0;

Xmatrix_alt = Xmatrix;
