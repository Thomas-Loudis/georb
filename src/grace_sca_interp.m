function [sca_int] = grace_sca_interp(acc,sca,scadpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grace_sca_interp:  GRACE star camera data interpolation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  GRACE star camera data interpolation to fit the rate of the
%  accelerometry data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - grace:    GRACE satellite ID
%             GRACE-A: grace=1, GRACE-B: grace=2 
% - acc:      accelerometry data in SRF (Science Reference Frame)
% - sca:      star camera assembly data (quatenions or Euler angles)
% - scadpint: Number of data points required for "sca" interpolation
%
% Output arguments:
% - facc:   Calibrated accelerometry array in inertial ICRS
%   facc = [MJD_TT tTT lin_accl_x lin_accl_y lin_accl_z]
%   MJDgps:      MJD in TT time scale including fraction of the day
%   tTT:         Seconds since 0h in TT time scale
%   lin_accl_i:  Calibrated linear acceleration along ICRS axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   April 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 22/04/2021, Thomas Papanikolaou
%             Code extracted from former grace_accproc.m function into an
%             individual function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial Accelerometry data array
% acc rate 1 sec
[sz1 sz2] = size(acc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Star camera interpolation
% sca rate 5 sec: interpolation at rate of 1 sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scadpint_prm = scadpint;
[sz5 sz6] = size(sca);
rate = 5;
sca_int = zeros(sz1,sz6);
k = 1;
jvaro = 0;
jvar = 0;
for i = 1 : sz1
    lp2 = 1;
    while lp2 == 1 && jvar <= sz5
        jvar = jvaro + 1;
        % Computing the Time difference because of different data rate
        if fix(acc(i,1)) - fix(sca(jvar,1)) == 0            
            if abs(acc(i,2) - sca(jvar,2)) < 10^-8
                sca_int(k,:) = sca(jvar,:);
                % next point
                k = k + 1;
                %break
                lp2 = 0;
                jvaro = jvar;
            elseif acc(i,2) - sca(jvar,2) < rate  % 5 sec
                % Lagrangian Interpolation
                [qo_int] = interp_Lagrange(sca(:,1),sca(:,3),acc(i,1),scadpint);
                [q1_int] = interp_Lagrange(sca(:,1),sca(:,4),acc(i,1),scadpint);
                [q2_int] = interp_Lagrange(sca(:,1),sca(:,5),acc(i,1),scadpint);
                [q3_int] = interp_Lagrange(sca(:,1),sca(:,6),acc(i,1),scadpint);           
                sca_int(k,:) = [acc(i,1) acc(i,2) qo_int q1_int q2_int q3_int];            
                clear qo_int q1_int q2_int q3_int
                % next point
                k = k + 1;
                %break
                lp2 = 0;
            end
        end
    end
end
% Same results are obtained if quaternion is multiplied by minus one (-1)
% sca_int = -1 * sca_int;
clear i jvar k rate lp2 jvaro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
