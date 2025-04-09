function [obs_reduced, obs_removed, Nepochs_remove, obs_residuals, rms_obs] = obs_pseudo_scan(sr1, sr2, sigma)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  obs_pseudo_scan.m
% 
% Purpose:
%  Pseudo_observations screening in order to remove observations with
%  residuals higher than the input threshold
% Statistics of the comparison between two Data Series Files.
%  The differences between the variables are computed in the common points
%  defined by the ID.
%  Statistical quantities are computed e.g. RMS, min, max..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - sr1: 1st Data Series file (used as reference data)
% - sr2: 2nd Data Series file
%
% Output arguments:
% - obs_reduced: Reduced array of reference observations (sr1) that pass the sigma criteria 
% - obs_removed: Removed observations from reference observations (sr1) that do not pass the sigma criteria 
% - RMS: Root Mean Square of each variable's differences
% - sr:  Array of sr1 & sr2 values at common points
%        sr = [id(common) sr1(common,2:end) sr2(common,2:end)]
% - min: Minimum value of each variable's differences
% - max: Maximum value of each variable's differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remarks:
% sr files format:     sr=[ID X1 X2 .... Xn]
% Differences are computed at common ID values:  dsr=[IDcommon dX1 ... dXn]
% RMS array format:  rms_dsr=[rms(dX1) ... rms(dXn)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Dr. Thomas Papanikolaou
% Created:  21 May 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Time_diglimit = 10^-8;

[sz1 sz2] = size(sr1);
[sz3 sz4] = size(sr2);
Ncol = sz2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute common epochs number 
epoch = 0;
wo = 1;
for i = 1 : sz3
    % 2nd Data series file
    id2 = sr2(i,1);
    % 1st (Reference) Data series file
    for w = wo : sz1
        id1 = sr1(w,1);
        if abs(id2 - id1) < Time_diglimit           
            epoch = epoch + 1; 
            wo = w + 1;
            break
        end
    end    
end
Nepochs = epoch;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dsr = zeros(Nepochs,Ncol);
% sr  = zeros(Nepochs, 1 + Ncol-1 + Ncol-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differences are computed at the common points defined by the ID
epoch = 0;
wo = 1;
Nepochs_remove = 0;
for i = 1 : sz3
    % 2nd Data series file
    id2 = sr2(i,1);
    % 1st (Reference) Data series file
    for w = wo : sz1
        id1 = sr1(w,1);
        if abs(id2 - id1) < Time_diglimit           
            epoch = epoch + 1;
            % Differences of each variable at common point
            dsr(epoch,:) = [sr2(i,:) - sr1(w,:)];
            dsr(epoch,1) = id1;
            sigma_test = 0;
            for j_delta = 2 : Ncol
                delta = sr2(i, j_delta) - sr1(w, j_delta);
                if abs(delta) >= sigma
                   sigma_test = 1; 
                end
            end
            if sigma_test == 1
                Nepochs_remove = Nepochs_remove + 1;
            end            
            % Array of sr1 and sr2 values at common points
            % sr(epoch,:) = [id1 sr1(w,2:end) sr2(i,2:end)];
            wo = w + 1;
            break
        end
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observations differences screening and remove epochs with differences
% higher than sigma threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Nepochs_remove == 0
   obs_reduced = sr1;
   obs_removed = 0;
   obs_residuals = dsr;
   Nepochs_reduced = Nepochs;
else
    
Nepochs_reduced = Nepochs - Nepochs_remove;
obs_reduced = zeros(Nepochs_reduced, Ncol);
obs_removed = zeros(Nepochs_remove, Ncol+Ncol-1);
obs_residuals = zeros(Nepochs_reduced, Ncol);

%[sz1 sz2] = size(sr1);
%[sz3 sz4] = size(sr2);
epoch = 0;
wo = 1;
i_epoch_removed = 0;
i_epoch_reduced = 0;
for i = 1 : sz3
    % 2nd Data series file
    id2 = sr2(i,1);
    % 1st (Reference) Data series file
    for w = wo : sz1
        id1 = sr1(w,1);
        if abs(id2 - id1) < Time_diglimit           
            epoch = epoch + 1;
            % Differences of each variable at common point
            delta_sr = [sr2(i,:) - sr1(w,:)];
            delta_residuals = delta_sr(2:end);
            sigma_test = 0;
            for j_delta = 2 : Ncol
                delta = sr2(i, j_delta) - sr1(w, j_delta);
                if abs(delta) >= sigma
                   sigma_test = 1; 
                end
            end
            if sigma_test == 1
                i_epoch_removed = i_epoch_removed + 1;
                obs_removed(i_epoch_removed,1:Ncol) = [sr1(w,:)];
                obs_removed(i_epoch_removed,1) = id1;                
                obs_removed(i_epoch_removed,Ncol+1:2*Ncol-1) = [delta_residuals(:)];
            end
            if sigma_test == 0
                i_epoch_reduced = i_epoch_reduced + 1;
                obs_reduced(i_epoch_reduced,:) = [sr1(w,:)];
                obs_reduced(i_epoch_reduced,1) = id1;                
                obs_residuals(i_epoch_reduced,:) = [sr2(i,:) - sr1(w,:)];
                obs_residuals(i_epoch_reduced,1) = id1;
            end            
            % Array of sr1 and sr2 values at common points
            % sr(epoch,:) = [id1 sr1(w,2:end) sr2(i,2:end)];
            wo = w + 1;
            break
        end
    end    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics
[sz5 sz6] = size(obs_residuals);
j = 1;
rms_obs = zeros(1,sz6-1);
for i = 2 : sz6
    rms_obs(1,j) = rms(obs_residuals(:,i));
    j = j + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
