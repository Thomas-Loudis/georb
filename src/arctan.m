function [angle] = arctan(y,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of an angle based on given tangent of the angle.
%
% angle = arctan(..)
% Checking the quadrant of angle that must be chosen.
% The computation is made counter-clockwise
%
% The result is computed in radians. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                      December 2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Basic angle a (refer to the first quadrant)
a = atan( abs( y/x ) );

% Reduction of the angle to the correct quadrant
if x > 0
    if y > 0
        angle = a;
    elseif y < 0
        angle = 2*pi-a;
    else
        angle = 0;
    end
elseif x < 0
    if y > 0
        angle = pi-a;
    elseif y < 0
        angle = pi+a;
    else
        angle = pi;
    end
else
    if y > 0
        angle = pi/2;
    elseif y < 0
        angle = 3*pi/2;
    end
end