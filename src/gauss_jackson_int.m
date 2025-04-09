function [Va] = gauss_jackson_int(zo,a,g,d,MSparam)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multistep numerical integration methods
% Gauss-Jackson or second sum methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Computation of initial values of first and second sums, V(^-1)ao and
%   V(^-2)ao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments
% - zo:  Initial State Vector in Celestial Reference System GCRS
%   zo = [zo1 zo2 zo3 zo4 zo5 zo6]'
%   zo1:  X component
%   zo2:  Y component
%   zo3:  Z component
%   zo4:  Vx component
%   zo5:  Vy component
%   zo6:  Vz component
% - a: Accelerations for first steps
%   ai = [aix aiy aiz]
% - MSparam:       Multistep method's parameters
%   MSparam(1,1):  ID number of selected multistep method
%   MSparam(2,1):  Order
%   MSparam(3,1):  Stepsize
%   MSparam(4,1):  Start integrator method for the first epochs
% - g: Coefficients gama g*j (gAM)
% - dc: Coefficients delta d*j (delta_c)
%      d*j = [d*0 d*1...d*j...d*m],  m: method's order 
%
% Output arguments
% - Va:  Total backwards differences including initial values of first and
%        second sums.
%        Va = [V2a
%              V1a
%              Va  ]
%        V2a = V(^-2)ao
%        V1a = V(^-1)ao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         May 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backwards Differences of acceleration ao
[N1, Na] = size(a);
%clear N1

% Va = zeros(Na,N1);
for i = 1 : Na
    [Vai] = backdiff(a(:,i));
    if i == 1
        [k1, k2] = size(Vai);
        Va = zeros(Na,k2);
    end
    Va(i,:) = Vai;
end
% clear i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial values of first and second sums
%  V(^-1)ao
sum_V1a = 0;
for j = 1 : MSparam(2,1)
    sum_V1a = sum_V1a + g(1,j+1) *  Va(:,j-1+1);
end
clear j
V1a = (1 / MSparam(3,1)) * [zo(4,1) zo(5,1) zo(6,1)]' - sum_V1a;
%  V(^-2)ao
sum_V2a = d(1,2) * V1a;
for j = 2 : MSparam(2,1)+1
    sum_V2a = sum_V2a + d(1,j+1) * Va(:,j-2+1);
end
clear j
V2a = (1 / MSparam(3,1)^2) * [zo(1,1) zo(2,1) zo(3,1)]' - sum_V2a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backwards Differences of acceleration ao
%  Total backwards differences including initial values for V(^-2)ao and
%  V(^-1)ao
Va = [V2a   V1a   Va];