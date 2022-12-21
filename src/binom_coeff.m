function [b] = binom_coeff(n,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binomial Coefficients
%
% Purpose:
%   Computation of the binomial coefficients that are reqired in
%   Adams-Bashforth integration method
%
% Input arguments
% - n,k:  Elements
%
% Output arguments
% - b:   Binomial coefficient of elements n,k according to the mathematical
%        scheme b_nk = (n k),
%        where n and k stand for the up and down elements respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         May 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


factor = 1;
for i = 0 : k-1
    factor = factor * (n - i);
end

b = factor / factorial(k);