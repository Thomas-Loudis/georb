function [yint] = interp_Lagrange(X,Y,xint,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lagrange Polynomials - Lagrangian Interpolation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Lagrangian interpolation is based on the computation of Lagrange
%  Polynomials and is suggested for data interpolation such as satellite
%  data, EOP (Earth Orientation Parameters), ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - X:     array of values of the independent variable
% - Y:     array of function values corresponding to X i.e. y=f(x)
% - xint:  X value for which interpolated estimate of y is required
% - dpint: number of data points used for interpolation
%
% Output arguments:
% - yout:  y value of the interpolation that refer to xint value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dpint_initial = dpint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test if interpolation is required or not
% Find index position of the first X value before xint 
[sz1 sz2] = size(X);
N_dataseries = sz1;

if dpint > N_dataseries
    fprintf('%s \n','Number of data points of interpolator is higher of the number of data series points. Interpolation number of data points is reduced accordingly.')
    dpint = N_dataseries;
end

dt_1epoch = X(sz1,1) - X(sz1-1,1);

if xint < X(1,1)
    fprintf('%s \n','Interpolation epoch is out of the data span')
    fprintf('%s %f   %s %f \n','Interpolation epoch :',xint,'Data series first epoch :',X(1,1))
    Xo_indx = 1 + (fix(dpint/2) - 1);
    interpolation = +1;
elseif xint > X(sz1,1) + dt_1epoch
    fprintf('%s \n','Interpolation epoch is out of the data span')
    fprintf('%s %f   %s %f \n','Interpolation epoch :',xint,'Data series last epoch :',X(sz1,1))        
    Xo_indx = sz1 -dpint+1 + (fix(dpint/2) - 1); 
    interpolation = +1;
end
for i = 1 : sz1
    if abs(xint - X(i,1)) < 10^-10
        yint = Y(i,1);
        interpolation = -1;
        break
    elseif xint - X(i,1) < 0
        Xo_indx = i-1;
        interpolation = +1;
        break
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if interpolation == 1
% Number of points before "xint"
Xno_before = Xo_indx;
% Number of points after "xint"
Xno_after = sz1 - Xo_indx;

% Define "dpint" limits before and after "xint" acording to the initial
% "dpint" ("dpint" is odd or even)
%dpint_initial = dpint;
dpint_limit1 = fix(dpint/2);
dpint_limit2 = dpint - fix(dpint/2);

% Special cases at the edges of data series (start and end)
% Start
if (dpint_limit1 > Xno_before) 
    Xo_indx = dpint_limit1;
end
% End 
if (dpint_limit2 > Xno_after)
    Xo_indx = sz1 - dpint_limit2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X,Y values of data points defined by the final value of "dpint"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find index position of the first X value for the data area defined
% by dpint (according if dpint is even or odd)
X1_indx = Xo_indx - (fix(dpint/2) - 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xdpint = zeros(dpint,1);
Ydpint = zeros(dpint,1);
for i = 1 : dpint
    Xdpint(i,1) = X(X1_indx+i-1,1);
    Ydpint(i,1) = Y(X1_indx+i-1,1);
end
% New X,Y matrices with the values for data points defined by dpint
X = Xdpint;
Y = Ydpint;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Langrange interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Data points
[n,m] = size(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of coefficients Li (xint)
% Prealocation
L = zeros(n,1);
for i = 1 : n
    L_numerator = 1 ;
    L_denominator = 1;
    for j = 1 : n
        if j == i
            % nothing
        else
            L_numerator = ( xint - X(j,1) ) * L_numerator ;
            L_denominator = ( X(i,1) - X(j,1) ) * L_denominator ;
        end
    end   
    L(i,1) = L_numerator / L_denominator ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of value of function Y(X) at point xint
% Estimation is realized by Interpolant Pn(xint) via Lagrange Polynomial
yint = 0;
for i = 1 : n
    yint = Y(i,1) * L(i,1) + yint;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
