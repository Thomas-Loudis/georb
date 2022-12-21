function [dp] = productdot(r1,r2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dot product between two vectors (scalar product)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - r1 : first vector   [x1 y1 z1]'
% - r2 : second vector  [x2 y2 z2]'
%
% Output arguments:
% - dp : Dot product of r1,r2 vectors (scalar)
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

dp = r1' * r2;
