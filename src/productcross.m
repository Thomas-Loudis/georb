function [cp] = productcross(r1,r2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross product between two vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - r1 : first vector   [x1 y1 z1]'
% - r2 : second vector  [x2 y2 z2]'
%
% Output arguments:
% - cp : Cross product of r1,r2 vectors  [x3 y3 z3]'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         June 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[n1 m1] = size(r1);
if n1 < m1
    r1 = r1';
end

[n2 m2] = size(r2);
if n2 < m2
    r2 = r2';
end

cp = [
r1(2,1)*r2(3,1) - r1(3,1)*r2(2,1)
r1(3,1)*r2(1,1) - r1(1,1)*r2(3,1)
r1(1,1)*r2(2,1) - r1(2,1)*r2(1,1)
];
