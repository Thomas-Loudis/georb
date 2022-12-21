function [RMS] = rms(dx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMS - root mean square
%
% Thomas D. Papanikolaou                                     December 2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[size_dx size2] = size(dx);
sum_dx2 = 0;
for i = 1 : size_dx
    sum_dx2 = sum_dx2 + dx(i,1)^2 ;
end

RMS = sqrt(sum_dx2 / size_dx) ;