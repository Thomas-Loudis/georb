function [accel_empirical] = empirical_cpr_init (EMP_FORCE_PARAM_yn)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: empirical_force_init.m (part of the old mainf_DOD.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Emprical Forces' Cycle-Per-Revolution coefficients initialisation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         June 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 20/04/2021, Thomas Papanikolaou
%             Code extracted and upgraded from the former mainf_DOD.m function
% 14/12/2022, Thomas Loudis Papanikolaou
%             Code upgrade to simplify the initialisation of the CPR
%             accelerations matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical forces matrix of parameters to be considered (yes/no)
[d1, d2] = size(EMP_FORCE_PARAM_yn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces parameters initialisation
% Preallocation to zero matrix
accel_empirical = zeros(d1,d2);

% Apriori values
emp_apriori = 10^-9;
accel_empirical = accel_empirical + emp_apriori;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Forces parameters: Number of unknown parameters to be estimated 
Nparam_EMP_FORCE = 0;
for i1 = 1 : d1
    for i2 = 1 : d2
        if EMP_FORCE_PARAM_yn(i1,i2) == 1
            Nparam_EMP_FORCE = Nparam_EMP_FORCE + 1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
