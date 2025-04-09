function [N_param_GRAV, Nparam_C, Nparam_S , C_degree_order, S_degree_order, dCnm_matrix, dSnm_matrix] = gravity_param_ic(degree_min, degree_max, order_min, order_max)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: prm_gfm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Gravity Field parameters to be estimated : Initialisation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Input confiugration file name *.in  
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                  15 April 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameter estimation
% grav_paramestim_degree_range(1,1) = degree_min;
% grav_paramestim_degree_range(1,2) = degree_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C coefficients
Ncoef = 0;
for n = degree_min : degree_max
    order_limit = order_max;
    if n < order_max
        order_limit = n;
    end
    for m = order_min : order_limit
        Ncoef = Ncoef + 1;
        C_degree_order(Ncoef,:) = [n m];
    end
end
Nparam_C = Ncoef;

% S coefficients
if order_min == 0
    % S_order_min = 1;
    S_order_min = 0;
end
if order_max == 0
    S_order_min = 0;
end
    S_order_min = 0;
Ncoef = 0;
for n = degree_min : degree_max
    order_limit = order_max;
    if n < order_max
        order_limit = n;
    end
    for m = S_order_min : order_limit
        Ncoef = Ncoef + 1;
        S_degree_order(Ncoef,:) = [n m];
    end
end
Nparam_S = Ncoef;

N_param_GRAV = Nparam_C + Nparam_S; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C coefficients matrix
dCnm_matrix = zeros(degree_max+1, degree_max+1);
dSnm_matrix = zeros(degree_max+1, degree_max+1);
init_value = 10^-15;
init_value = 0;
% Cnm estimated corrections
for i_C = 1 : Nparam_C
    n = C_degree_order(i_C,1);
    m = C_degree_order(i_C,2);
    dCnm_matrix(n+1,m+1) = init_value;
end
% Snm estimated corrections
for i_C = 1 : Nparam_S
    n = S_degree_order(i_C,1);
    m = S_degree_order(i_C,2);
    if m == 0
    dSnm_matrix(n+1,m+1) = 0;
    else 
    dSnm_matrix(n+1,m+1) = init_value;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
