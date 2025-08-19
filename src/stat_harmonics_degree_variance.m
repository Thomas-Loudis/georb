function [degree_variance, degree_amplitude, degree_variance_sqrt, degree_amplitude_sqrt] = stat_harmonics_degree_variance (C_matrix, S_matrix, n_max, n_min)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: stat_harmonics_degree_variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Statistical assessment of estimated gravity field models
%  Compute degree amplitude and degree variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                   28 June 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n_trunc = n_max;
% n_min = 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics :: Degree Variances and Amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
degree_variance  = zeros(n_trunc+1,1);
degree_amplitude = zeros(n_trunc+1,1);
for n = n_min : n_trunc
    sigma2_sum = 0;
    for m = 0 : n
        sigma2_sum = sigma2_sum + (C_matrix(n+1,m+1)^2 + S_matrix(n+1,m+1)^2);
    end
    degree_variance(n+1,1)  = sigma2_sum;
    degree_amplitude(n+1,1) = sigma2_sum / (2*n+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

degree_variance_sqrt = sqrt(degree_variance);

degree_amplitude_sqrt = sqrt(degree_amplitude);
