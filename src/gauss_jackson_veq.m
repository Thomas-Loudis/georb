function [veqZ] = gauss_jackson_veq(V_veqZr_dv2,g,d,MSparam)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gauss-Jackson (second sum) multistep numerical integration method
% for integration of Variational Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Numerical integration of second-order differential equation.
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
% - Va: Total backwards differences including initial values of first and
%       second sums.
%       Va = [V2a
%             V1a
%             Va  ]
%       V2a = V(^-2)ao
%       V1a = V(^-1)ao
% - MSparam:       Multistep method's parameters
%   MSparam(1,1):  ID number of selected multistep method
%   MSparam(2,1):  Order
%   MSparam(3,1):  Stepsize
%   MSparam(4,1):  Start integrator method for the first epochs
% - g: Coefficients gama gj or g*j
% - d: Coefficients delta dj (delta_s) or d*j (delta_c)
%      dj = [d0 d1...dj...dm],  m: method's order 
%
% Output arguments
% - z:  Solution at epoch t=to+h in GCRS
%   z = [z1 z2 z3 z4 z5 z6]'
%   z1:  X component
%   z2:  Y component
%   z3:  Z component
%   z4:  Vx component
%   z5:  Vy component
%   z6:  Vz component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         May 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution of VEQ Integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% veqZr array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz1 sz2 sz3] = size(V_veqZr_dv2);
for i1 = 1 : sz1
    for i2 = 1 : sz2
        sumV_veqZr_dv2 = 0;
        for j = 0 : MSparam(2,1)+1
            sumV_veqZr_dv2 = sumV_veqZr_dv2 + d(1,j+1) * V_veqZr_dv2(i1,i2,j-2+3); 
        end
        sumV_array(i1,i2) = sumV_veqZr_dv2;
    end
end
veqZr = MSparam(3,1)^2 * sumV_array;
clear sz1 sz2 sz3 i1 i2 i3 sumV_veqZr_dv2 sumV_array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% veqZv array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz1 sz2 sz3] = size(V_veqZr_dv2);
for i1 = 1 : sz1
    for i2 = 1 : sz2
        sumV_veqZr_dv2 = 0;
        for j = 0 : MSparam(2,1)
            sumV_veqZr_dv2 = sumV_veqZr_dv2 + g(1,j+1) * V_veqZr_dv2(i1,i2,j-1+3);
        end
        sumV_array(i1,i2) = sumV_veqZr_dv2;
    end
end
veqZv = MSparam(3,1) * sumV_array;
clear sz1 sz2 sz3 i1 i2 i3 sumV_veqZr_dv2 sumV_array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% veqZ array 
veqZ = [veqZr
        veqZv];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 