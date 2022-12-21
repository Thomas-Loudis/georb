function [orb2] = mjd2mjdtt(orb,fixmjd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJD is including the fraction of the day
%
% MJD: convert or not to integer
% fixmjd = 1 : convert to integer
% fixmjd = 0 : do not convert to integer
%
% Time conversion: MJD column conversion to 2 columns, MJD and t (sec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou,  AUTH                                    May 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[sz1, sz2] = size(orb);
orb2 = zeros(sz1,sz2+1);
for i = 1 : sz1
    mjdi = orb(i,1);
    [ti,D,M,Y] = MJD_inv(mjdi);
    if fixmjd == 1
        mjdi = fix(mjdi);
    end
    orb2(i,:) = [mjdi ti orb(i,2:end)];
end
