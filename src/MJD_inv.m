function [sec,day,month,year] = MJD_inv(MJD)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion between the Modified Julian Day Number (MJD) and the Calendar
% Date
%
% Conversion from Modified Julian day number to Civil date.
% Civil Date is expressed in Gregorian calendar.
%
% The algorithm is based on the algorithms given by Fliegel and van
% Flandern (1968).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - MJD: Modified julian day number
%
% Output arguments:
% - sec: Seconds since 00h, fraction of the day
% - D: number of day
% - M: number of month
% - Y: Year
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference: 
% Fliegel, H. F. & van Flandern, T. C. 1968,  A Machine Algorithm for
% Processing Calendar Dates, Communications of the ACM, 11:657.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remarks:
%  The function entier(x) or briefly [x] is defined as as the smallest
%  integer which is smaller than or equal to x, i.e. [x] <= x < x+1.
%  For positive numbers it is equal to the integral part of x. For negative
%  (non-integer) numbers it is equal to the integral part of x minus one.
%  Here, Matlab's function floor(x) is used for this purpose.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                   November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 30/10/2022, Dr. Thomas Papanikolaou
%             Algorithm has been added 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Julian Day number 
JD = MJD + 2400000.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alg = 2;


if alg == 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of civil date by intermediate auxiliary quantities
% a,b,c,d,e,f
q = MJD - floor(MJD);

a = floor(MJD) + 2400001;

if a<2299161
    b = 0;
else
    b = floor( (a-1867216.25) / 36524.25 );
end

if a<2299161
    c = a + 1524;
else
    c = a + b - floor(b/4) + 1525;
end

%d = floor( (c-121.1) / 365.25);
d = floor( (c-122.1) / 365.25);
e = floor(365.25 * d);
f = floor( (c-e) / 30.6001);

% compute fraction of the day - UT in seconds
sec = q * 24 * 3600;

% compute Day (D)
day = c - e - floor(30.6001 * f) ;
%day = c - e - floor(30.6001*f) + q

% compute Month (M)
month = f - 1 - 12 * floor(f/14);

% compute Year(Y)
year = d - 4715 - floor((7+month)/10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif alg == 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = fix(JD + 0.5) + 1537;
c = fix( (b - 122.1) / 365.25);
d = fix( 365.25 * c);
e = fix( (b-d) / 30.6001);

day   = b - d - fix(30.6001 * e);
month = e - 1 - 12 * fix(e / 14);
year  = c - 4715 - fix( (7 + month) / 10);

sec = (MJD - fix(MJD)) * 86400;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif alg == 3
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Julian Day number 
JD_full = JD;
JD_day = fix(JD);
JD = JD_day;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fliegel and van Flandern (1968)
L  = JD + 68569;
N  = (4 * L) / 146097;
L  = L - (146097 * N + 3) / 4;
I  = (4000 * (L+1) ) / 1461001;
L  = L - (1461 * I) / 4 + 31;
J  = (80 * L) / 2447;
K  = L - (2447 * J) / 80;
L  = J / 11;
J  = J + 2 - 12 * L;
I  = 100 * (N - 49) + I + L;

day   = K;
month = J;
year  = I;

Sec_00 = ( MJD - fix(MJD) ) * 86400;
sec = Sec_00;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


elseif alg == 4

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Julian Day number 
JD = MJD + 2400000.5;
JD_full = JD;
JD_day = fix(JD);
JD = JD_day;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = JD + 68569;
N = ( 4*L ) / 146097;
L = L - ( 146097*N + 3 ) / 4;
I = ( 4000 * (L+1) ) / 1461001;
L = L - ( 1461*I ) / 4 + 31;
K = ( 80*L ) / 2447;
ID = L - ( 2447*K ) / 80;
L = K / 11;
IM = K + 2 - 12*L;
IY = 100 * ( N-49 ) + I + L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

day   = ID;
month = IM;
year  = IY;

Sec_00 = ( MJD - fix(MJD) ) * 86400;
sec = Sec_00;

end     
      
      