function [JD,MJD] = MJD_date(t,D,M,Y)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion between the Modified Julian Day Number (MJD) and the Calendar
% Date.
%
% Conversion from Civil date to Julian day number and MJD.
% Civil Date is expressed in Gregorian calendar.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - t: time in seconds since 0 hours
% - D: number of day
% - M: number of month
% - Y: Year
%
% Output arguments:
% - JD:  Julian day number
% - MJD: Modified Julian day number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Day including the fraction of the day
Dfr = D + t / (24*3600);

if M <= 2
    y = Y - 1;
else
    y = Y;
end

if M <= 2
    m = M + 12;
else
    m = M;
end

if Y < 1582
    B = -2 + floor((Y + 4716) / 4) - 1179;
elseif Y > 1582
    B = floor(Y/400) - floor(Y/100) + floor(Y/4);
    %B = -15 + floor(Y/4);
elseif Y == 1582
    if M < 10
        B = -2 + floor((Y + 4716) / 4) - 1179;
    elseif M > 10
        B = floor(Y/400) - floor(Y/100) + floor(Y/4);
    elseif M == 10
        if D <= 4
            B = -2 + floor((Y + 4716) / 4) - 1179;
        elseif D >= 10
            B = floor(Y/400) - floor(Y/100) + floor(Y/4);
        elseif D > 4 && D < 10
            %fprintf('%s \n', 'Date is between 4 and 10 Ocotber 1582')
            B = 0;
            [JDt,MJDt] = mjd3(t,D,M,Y);
        end
    end
end

if B == 0
    % Alternative 2nd algorthm (function mjd3.m)
    MJD = MJDt;
    JD = JDt;
else
    MJD = 365*y -679004 + B + floor(30.6001*(m+1)) + Dfr;
    %MJDpt1 = 365*y -679004 + B
    %MJDpt2 = floor(30.6001*(m+1))
    JD = MJD + 2400000.5;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
