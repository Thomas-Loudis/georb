function [y,eABm] = AB(yo,f,MSparam,g,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multistep numerical integration methods
% Adams-Bashforth method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Numerical integration of first-order differential equation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments
% - yo:  Initial State Vector in Celestial Reference System GCRS
%   yo = [yo1 yo2 yo3 yo4 yo5 yo6]'
%   yo1:  X component
%   yo2:  Y component
%   yo3:  Z component
%   yo4:  Vx component
%   yo5:  Vy component
%   yo6:  Vz component
% - f: Fuction evaluations
%   fi = [fi_vx fi_vy fi_vz fi_ax fi_ay fi_az]
% - MSparam:       Multistep method's parameters
%   MSparam(1,1):  ID number of selected multistep method
%   MSparam(2,1):  Order
%   MSparam(3,1):  Stepsize
%   MSparam(4,1):  Start integrator method for the first epochs
% - g: Coefficients gama gj
% - b: Coefficients bmj
%      bmj = [bm1...bmj...bmm],  m: method's order 
%
% Output arguments
% - y:  Solution at epoch t=to+h in GCRS
%   y = [y1 y2 y3 y4 y5 y6]'
%   y1:  X component
%   y2:  Y component
%   y3:  Z component
%   yo:  Vx component
%   yo:  Vy component
%   yo:  Vz component
% - eABm:  local truncation error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         May 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local variables for the multistep method's parameters :: Order 
ord = MSparam(2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Increment function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summarized Formulas (Substituting the backwards differences)
FAB = 0;
for j = 1 : ord
    FAB = FAB + b(1,j) * f(j,:)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backwards differences
[n1, Nf] = size(f);
clear n1
for i = 1 : Nf
    [Vfi] = backdiff(f(:,i));
    Vf(i,:) = Vfi;
end
clear i
FAB_bw = 0;
for j = 0 : ord-1
    FAB_bw = FAB_bw + g(1,j+1) * Vf(:,j+1);
end
clear j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution of integration
% 1. Summarized Formulas
% 2. Backwards differences
%FAB = FAB_bw;
h = MSparam(3,1);
y = yo + h * FAB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local truncation error
%eABm = h * abs( g(1,ord+1) * Vf(1,ord+1) )  % Vmf of order m is required 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%