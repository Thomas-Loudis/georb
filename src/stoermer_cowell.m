function [r] = stoermer_cowell(ri,a,MSparam,d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multistep numerical integration methods
% Stoermer and Cowell methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Numerical integration of second-order differential equation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments
% - ri:  Initial Position Vector in Celestial Reference System GCRS
%        for two points at initial epoch 0 (i) and epoch -1 (i-1).
%   ri = [X(i-1) Y(i-1) Z(i-1)
%         X(i)   Y(i)   Z(i)  ]
% - a: Accelerations for first steps
%   ai = [aix aiy aiz]
% - MSparam:       Multistep method's parameters
%   MSparam(1,1):  ID number of selected multistep method
%   MSparam(2,1):  Order
%   MSparam(3,1):  Stepsize
%   MSparam(4,1):  Start integrator method for the first epochs
% - d: Coefficients delta dj (delta_s) or d*j (delta_c)
%      dj = [d0 d1...dj...dm],  m: method's order 
%
% Output arguments
% - r:  Solution at epoch t=to+h in GCRS
%   r = [y1 y2 y3]'
%   y1:  X component
%   y2:  Y component
%   y3:  Z component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         May 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backwards differences of accelaration ai
[n1 Na] = size(a);
clear n1
for i = 1 : Na
    [Vai] = backdiff(a(:,i));
    Va(i,:) = Vai;
end
clear i Na
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution of integration
sumVa = 0;
for j = 0 : MSparam(2,1)-1
    sumVa = sumVa +  d(1,j+1) * Va(:,j+1);
end
clear j
h = MSparam(3,1);
r = 2 * ri(2,:)' - ri(1,:)' + h^2 * sumVa;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%