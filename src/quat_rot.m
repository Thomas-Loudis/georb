function [R] = quat_rot(q)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quat_rot:  Quaternions Rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the rotation matrix based on quaternion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - q:   quaternion
%
% Output arguments:
% - R:   Rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
% Quaternion form of input argument "q":
% q = [qo q1 q2 q3]'
% q = qo + q1*i + q2 * j + q3 * k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                 October 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quaternion form: q = [qo q1 q2 q3]'
qo = q(1,1);
q1 = q(2,1);
q2 = q(3,1);
q3 = q(4,1);

% Rotation matrix based on quaternion
R = [ qo^2+q1^2-q2^2-q3^2	2*(q1*q2+qo*q3)		2*(q1*q3-qo*q2)
	2*(q1*q2-qo*q3)		qo^2-q1^2+q2^2-q3^2	2*(q2*q3+qo*q1)
	2*(q1*q3+qo*q2)		2*(q2*q3-qo*q1)		qo^2-q1^2-q2^2+q3^2 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
