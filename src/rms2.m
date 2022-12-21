function [RMS] = rms(dx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMS - root mean square
%
% Thomas D. Papanikolaou                                     December 2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 13/08/2012   Computations are now indepedent of dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[sz1 sz2] = size(dx);
if sz1 >= sz2
    size_dx = sz1;
    dm = 1;
else
    size_dx = sz2;
    dm = 2;
end
   

sum_dx2 = 0;
for irun = 1 : size_dx
    if dm == 1
        sum_dx2 = sum_dx2 + dx(irun,1)^2 ;
    elseif dm == 2
        sum_dx2 = sum_dx2 + dx(1,irun)^2 ;
    end
end

RMS = sqrt(sum_dx2 / size_dx) ;