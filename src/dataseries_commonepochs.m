function [epochs_common, sr1_common, sr2_common] = dataseries_commonepochs(sr1,sr2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  compstat.m
% 
% Purpose:
%  Statistics of the comparison between two Data Series Files.
%  The differences between the variables are computed in the common points
%  defined by the ID.
%  Statistical quantities are computed e.g. RMS, min, max..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - sr1: 1st Data Series file (used as reference data)
% - sr2: 2nd Data Series file
%
% Output arguments:
% - delta_series: Differences between the input files
%   >> delta_series = [id(common) sr2-sr1(:,2:end)]
% - rms_dsr: Root Mean Square of the individual variable's differences
%   rms_dsr = RMS(delta_series)
% - epochs_common:  Array of common epochs between sr1 & sr2 series 
%   series_common = [id(common) sr1(common,2:end) sr2(common,2:end)]
% - min: Minimum value of each variable's differences
% - max: Maximum value of each variable's differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remarks:
% sr files format:     sr=[ID X1 X2 .... Xn]
% Differences are computed at common ID values:  dsr=[IDcommon dX1 ... dXn]
% RMS array format:  rms_dsr=[rms(dX1) ... rms(dXn)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                           May 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Time_diglimit = 10^-8;
[sz1 sz2] = size(sr1);
[sz3 sz4] = size(sr2);
Ncol = sz2;

% Matrices Initialisation / Preallocation
Nepochs_min = sz1;
if sz3 < sz1
    Nepochs_min = sz3;
end
dsr = zeros(Nepochs_min,1);

sr1_common_temp = zeros(Nepochs_min,sz2);
sr2_common_temp = zeros(Nepochs_min,sz4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differences are computed at the common points defined by the ID
epoch = 0;
wo = 1;
for i = 1 : sz3
    % 2nd Data series file
    id2 = sr2(i,1);
    % 1st (Reference) Data series file
    for w = wo : sz1
        id1 = sr1(w,1);
        % if abs(id2 - id1) < Time_diglimit
        delta_id = abs(id2 - id1);
        if delta_id < Time_diglimit
            epoch = epoch + 1;
            % Differences of each variable at common point
            dsr(epoch,1) = id1;
            % Array of sr1 and sr2 values at common points
            sr1_common_temp(epoch,:) = sr1(w,:);
            sr2_common_temp(epoch,:) = sr2(i,:);
             
            wo = w + 1;
            break
        end
    end    
end
Nepochs_common = epoch;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear empty data entries
epochs_common = dsr(1:Nepochs_common,:);
sr1_common    = sr1_common_temp(1:Nepochs_common,:);
sr2_common    = sr2_common_temp(1:Nepochs_common,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

