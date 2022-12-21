function [Vall_veqZr_dv2] = gauss_jackson_int_veq(veqZo,veqZr_dv2,g,d,MSparam)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gauss-Jackson (second sum) multistep numerical integration method
% for integration of Variational Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%   Computation of initial values of first and second sums, V(^-1)ao and
%   V(^-2)ao required by Gauss-Jackson integrator formula
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
% Thomas D. Papanikolaou                                         May 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial values veqZr, veqZv
veqZr = veqZo(1:3,:);
veqZv = veqZo(4:6,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backwards Differences of veqZr_dv2 (3D array)
[sz1 sz2 sz3] = size(veqZr_dv2);
for i1 = 1 : sz1
    for i2 = 1 : sz2
        jrun = 0;
        for i3 = 1 : sz3
            jrun = jrun + 1;
            veqZr_dv2_i1i2(jrun,1) = veqZr_dv2(i1,i2,i3);
        end
        clear i3 jrun
        % veqZr_dv2_i1i2   : column matrix
        % V_veqZr_dv2_i1i2 : row matrix
        [V_veqZr_dv2_i1i2] = backdiff(veqZr_dv2_i1i2);
        for i3 = 1 : sz3
            V_veqZr_dv2(i1,i2,i3) = V_veqZr_dv2_i1i2(1,i3);
        end       
        clear i3 V_veqZr_dv2_i1i2 veqZr_dv2_i1i2
    end
end
clear i1 i2 i3 sz1 sz2 sz3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial values of first and second sums
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  V(^-1)ao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz1 sz2 sz3] = size(V_veqZr_dv2);
for i1 = 1 : sz1
    for i2 = 1 : sz2
        sum_V1 = 0;
        for j = 1 : MSparam(2,1)
            sum_V1 = sum_V1 + g(1,j+1) *  V_veqZr_dv2(i1,i2,j-1+1);
        end
        sum_V1_veqZr_dv2(i1,i2) = sum_V1;
    end 
end
clear j sz1 sz2 sz3
V1_veqZr_dv2 = (1 / MSparam(3,1)) * veqZv - sum_V1_veqZr_dv2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  V(^-2)ao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_V2_d1 = d(1,1+1) * V1_veqZr_dv2;
[sz1 sz2 sz3] = size(V_veqZr_dv2);
for i1 = 1 : sz1
    for i2 = 1 : sz2
        sum_V2 = 0;
        for j = 2 : MSparam(2,1)+1
            sum_V2 = sum_V2 + d(1,j+1) * V_veqZr_dv2(i1,i2,j-2+1);  
        end
        sum_V2_veqZr_dv2(i1,i2) = sum_V2 + sum_V2_d1(i1,i2);
    end
end
clear j sz1 sz2 sz3
V2_veqZr_dv2 = (1 / MSparam(3,1)^2) * veqZr - sum_V2_veqZr_dv2;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backwards Differences of acceleration ao
%  Total backwards differences including initial values for V(^-2)ao and
%  V(^-1)ao
[sz1 sz2 sz3] = size(V_veqZr_dv2);
Vall_veqZr_dv2 = zeros(sz1,sz2,sz3+2);
Vall_veqZr_dv2(:,:,1)         = V2_veqZr_dv2;
Vall_veqZr_dv2(:,:,2)         = V1_veqZr_dv2;
Vall_veqZr_dv2(:,:,3 : sz3+2) = V_veqZr_dv2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
