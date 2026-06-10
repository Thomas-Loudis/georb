function [Amatrix] = neq_veq2Amatrix (veqZarray,veqParray, epochs_common, epochs_veq) 

%-------------------------------------------------------------------------- 
% Copyright (C) 2007-present  Thomas (Loudis) Papanikolaou  
%  
% GEORB :: 'Gravity and Precise Orbit Determination system' 
%  
% GEORB is a software for Precise Orbit Determination (POD), gravity field  
% recovery and mission design 
%  
% GEORB is licensed under the GNU Affero General Public License v3.0 while  
% for information regarding dual-license options visit the official website 
% www.georb.gr  
%-------------------------------------------------------------------------- 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Design matrix formulation based on the partial derivatives of the 
% Variational Equations' solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input arguments: 
% - veqZarray:      Variational Equations (VEQ) solution for the initial state 
%                   vector (State transition matrix) per epoch                    
% - veqParray:      Variational Equations (VEQ) solution for additional parameters 
%                   (force related parameters) per epoch  
% - epochs_common:  Common epochs matrix in MJD including fraction of the day 
% - epochs_veq:     VEQ matrices' Epochs matrix in MJD including fraction of the day 
% 
% Output arguments: 
% - Amatrix:        Design matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dr. Thomas Loudis Papanikolaou                             11 August 2025 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
[n1_epochs_veq,n2_epochs_veq] = size(epochs_veq); 
[N_commonepochs, n2] = size(epochs_common); 
 
[n1_VEQ_Z,n2_VEQ_Z] = size(veqZarray); 
[n1_VEQ_P,n2_VEQ_P] = size(veqParray); 
 
% obstype  :  1/2  [r/v] 
obstype = 1; 
if obstype == 1 
    Nobsset = 3; 
elseif obstype == 2 
    Nobsset = 6; 
end 
 
% Matrices Preallocation 
Amatrix_Z = zeros(N_commonepochs * Nobsset, n2_VEQ_Z - 1); 
Amatrix_P = zeros(N_commonepochs * Nobsset, n2_VEQ_P - 1); 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Design Matrix :: formulated at the common epochs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Nepochs = 0; 
i0_ref = 1; 
for i_epochs = 1 : N_commonepochs 
    ti_common = epochs_common(i_epochs,1); 
    for iref = i0_ref : n1_epochs_veq 
        % Reference Orbit 
        ti_ref = epochs_veq(iref,1); 
        if abs(ti_ref - ti_common) < 10^-8 
            Nepochs = Nepochs + 1; 
            % timearray = [ti_ref; ti_ref; ti_ref];                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % Design Matrix values obtained by VEQ arrays values 
            % veqZarray 
            Amatrix_ti = veqZarray( (iref-1)*6+1 : (iref-1)*6+Nobsset , 2:end); 
            Amatrix_Z((Nepochs-1)*Nobsset+1 : Nepochs*Nobsset,:) = Amatrix_ti; 
            % Amatrix_time((Nepochs-1)*Nobsset+1 : Nepochs*Nobsset,:) = [timearray Amatrix_ti]; 
 
            % veqParray 
            if size(veqParray,1) > 1  
            AmatrixP_ti = veqParray( (iref-1)*6+1 : (iref-1)*6+Nobsset , 2:end); 
            Amatrix_P((Nepochs-1)*Nobsset+1 : Nepochs*Nobsset,:) = AmatrixP_ti; 
            end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            i0_ref = iref; 
        end 
    end 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Design Matrix: Final Amatrix  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if size(veqParray,1) > 1  
Amatrix_ZP = [Amatrix_Z Amatrix_P]; 
else 
Amatrix_ZP = Amatrix_Z; 
end 
 
Amatrix = Amatrix_ZP; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
