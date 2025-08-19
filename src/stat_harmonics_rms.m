function [rms_CS] = stat_harmonics_rms (CS_matrix, n_min, n_max, m_min)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: stat_harmonics_rms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Compute RMS of spherical harmonic coefficients numerical differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                   28 June 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics :: RMS of delta coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_CS = CS_matrix;
N_delta = ( (n_max+1)^2 - (n_max+1) ) / 2 + (n_max+1);
delta_matrix = zeros(N_delta,1);

i_delta = 0;
for n = n_min : n_max
    for m = m_min : n
        i_delta = i_delta + 1;
        delta_matrix(i_delta,1) = delta_CS(n+1,m+1);
    end
end
% i_delta
delta_matrix_out = delta_matrix(1:i_delta,1);

% RMS of delta coefficients
rms_CS = rms(delta_matrix_out);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%