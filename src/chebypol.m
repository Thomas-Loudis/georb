function [ft,dft] = chebypol(coeff,t,t1,t2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chebychev polynomials approximation
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Evaluation of a function's value and its derivative based on Chebychev
%  approximation via known Chebychev coefficients.
%
% Input arguments
% - coeff : Chebychev coefficients matrix (column matrix)   
% - t     : Computation (interpolation) time
% - t1    : Start time of Chebychev coefficients
% - t2    : End time of Chebychev coefficients
%
%   t,t1,t2 are given in TT (Terrestrial Time) in the form of Julian Date
%   number (JD).
%
% Output arguments:
% - ft   : Function at interpolation epoch t in ...  (e.g. Planet position)
% - dft  : Derivative of function ft at t in ...     (e.g. Planet velocity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remarks
% Equations index "i" varies from 0 to n-1 while the respective index for
% the matrix elements varies from 1 to n. 
% Example:
%  Chebyshev coefficients ao,a1...an-1
%  Matrix of Chebyshev coefficients coeff(1,1),coeff(2,1).....coeff(n,1) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                               November  2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 16/03/2025    Dr. Thomas Loudis Papanikolaou
%               Code upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coeff_no = length(coeff);
% Polynomial order (o ... n-1)
n = coeff_no;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% taph (t in greek) computation (t mapped to interval [-1 1])
taph = 2 * (t - t1) / (t2 - t1) - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chebychev polynomials T : Computation via recurence formula
To = 1;
T1 = taph;
T(0 +1,1) = To;
T(1 +1,1) = T1;
for i = 1 : n
   T(i+1 +1,1) = 2 * taph * T(i +1,1) - T(i-1 +1,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chebychev approximation evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation based on 2 ways:
% 1. Summ equation
% 2. Clenshaw algorithm
% 
% Select by setting variable "cheby_eval" with the respective ID
cheby_eval = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cheby_eval == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Summ equation
ft = 0;
for i = 0 : n-1
    ft = ft + coeff(i +1,1) * T(i +1,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif cheby_eval == 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 2. Clenshaw algorithm
f(n +1,1) = 0;
f(n+1 +1,1) = 0;
for i = n-1 : -1 : 1
    f(i +1,1) = 2 * taph * f(i+1 +1,1) - f(i+2 +1,1) + coeff(i +1);
end
ft = taph * f(1 +1,1) - f(2 +1,1) + coeff(0 +1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative of Chebychev approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivatives of Chebychev coefficients
coeff_dev(n-1 +1,1) = 0;
coeff_dev(n +1,1) = 0;
for i = n-2 :-1: 1
    coeff_dev(i +1,1) = coeff_dev(i+2 +1,1) + 2 * (i+1) * coeff(i+1 +1);
end
coeff_dev(0 +1,1) = coeff(1 +1,1) + coeff_dev(2 +1,1) / 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative of the function evaluation via Chebychev approximation
dft = 0;
for i = 0 : n-2
    dft = dft + coeff_dev(i +1,1) * T(i +1,1);
end
dft = (2 / (t2-t1) ) * dft;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
