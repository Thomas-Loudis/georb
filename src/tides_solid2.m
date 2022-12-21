function [dCnm,dSnm] = tides_solid2(mjd,eop,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tides_solid1 : Solid Earth Tides
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Solid Earth Tides based on IERS Conventions 2010. This function impleme-
%  nts the Step2 according to the IERS Conventions 2010 which refer to the 
%  computation of the frequency depedent corrections to the  geopotential's
%  spherical harmonic coefficients.
%
% Input arguments
% - mjd     : MJD in Terrestrial Time (TT) scale including the fraction of
%             the day 
% - eop:    Earth Orientation Parameters (EOP) data that are required for
%           the orbit arc length
% - dpint:  Number of data points (days) that are required for the EOP
%           interpolation to the computation epoch
%
% Output arguments:
% - dCnm    : Cnm corrections matrix
% - dSnm    : Snm corrections matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
%  Computed dCnm and dSnm are stored into lower triangular matrices.
%  Coefficient dCnm corresponds to matrix element dCnm(n+1,m+1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                   June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% dCnm = zeros(5,5);
% dSnm = zeros(5,5);

% Delaunay variables (in radians)
[F1,F2,F3,F4,F5] = delaunay_variables(mjd);

% Greenwich Mean Sidereal Time (GMST) in radians
[thetag] = iers_gmst(mjd,eop,dpint);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dC20 correction
% Table 6.5b
% Matrix format : Columns : Doodson_No. l l' F D Omega Amp.(ip) Amp.(op)
table65b  = [
55.565  0  0  0  0  1  16.6  -6.7
55.575  0  0  0  0  2  -0.1   0.1
56.554  0 -1  0  0  0  -1.2   0.8
57.555  0  0 -2  2 -2  -5.5   4.3
57.565  0  0 -2  2 -1   0.1  -0.1
58.554  0 -1 -2  2 -2  -0.3   0.2
63.655  1  0  0 -2  0  -0.3   0.7
65.445 -1  0  0  0 -1   0.1  -0.2
65.455 -1  0  0  0  0  -1.2   3.7
65.465 -1  0  0  0  1   0.1  -0.2
65.655  1  0 -2  0 -2   0.1  -0.2
73.555  0  0  0 -2  0   0.0   0.6
75.355 -2  0  0  0  0   0.0   0.3
75.555  0  0 -2  0 -2   0.6   6.3
75.565  0  0 -2  0 -1   0.2   2.6
75.575  0  0 -2  0  0   0.0   0.2
83.655  1  0 -2 -2 -2   0.1   0.2
85.455 -1  0 -2  0 -2   0.4   1.1
85.465 -1  0 -2  0 -1   0.2   0.5
93.555  0  0 -2 -2 -2   0.1   0.2
95.355 -2  0 -2  0 -2   0.1   0.1
 ];         
% dC20 computation
[sz1 sz2] = size(table65b);
dC20 = 0;
for i = 1 : sz1
    % thetaf (in radians)
    thetaf = - table65b(i,2:6) * [F1 F2 F3 F4 F5]';
    dC20 = dC20 + table65b(i,7) * 10^(-12) * cos(thetaf) -  table65b(i,8) * 10^(-12) * sin(thetaf);
    clear thetaf
end
clear sz1 sz2 i thetaf
dCnm(2+1,0+1) = dC20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dC21,dS21 corrections 
% Table 6.5a
% Matrix format : Columns : Doodson_No. l l' F D Omega Amp.(ip) Amp.(op)
table65a  = [
125.755  2  0  2  0  2   -0.1    0.0
127.555  0  0  2  2  2   -0.1    0.0
135.645  1  0  2  0  1   -0.1    0.0
135.655  1  0  2  0  2   -0.7    0.1
137.455 -1  0  2  2  2   -0.1    0.0
145.545  0  0  2  0  1   -1.3    0.1
145.555  0  0  2  0  2   -6.8    0.6
147.555  0  0  0  2  0    0.1    0.0
153.655  1  0  2 -2  2    0.1    0.0
155.445 -1  0  2  0  1    0.1    0.0
155.455 -1  0  2  0  2    0.4    0.0
155.655  1  0  0  0  0    1.3   -0.1
155.665  1  0  0  0  1    0.3    0.0
157.455 -1  0  0  2  0    0.3    0.0
157.465 -1  0  0  2  1    0.1    0.0
162.556  0  1  2 -2  2   -1.9    0.1
163.545  0  0  2 -2  1    0.5    0.0
163.555  0  0  2 -2  2  -43.4    2.9
164.554  0 -1  2 -2  2    0.6    0.0
164.556  0  1  0  0  0    1.6   -0.1
165.345 -2  0  2  0  1    0.1    0.0
165.535  0  0  0  0 -2    0.1    0.0
165.545  0  0  0  0 -1   -8.8    0.5
165.555  0  0  0  0  0  470.9  -30.2
165.565  0  0  0  0  1   68.1   -4.6
165.575  0  0  0  0  2   -1.6    0.1
166.455 -1  0  0  1  0    0.1    0.0
166.544  0 -1  0  0 -1   -0.1    0.0
166.554  0 -1  0  0  0  -20.6   -0.3
166.556  0  1 -2  2 -2    0.3    0.0
166.564  0 -1  0  0  1   -0.3    0.0
167.355 -2  0  0  2  0   -0.2    0.0
167.365 -2  0  0  2  1   -0.1    0.0
167.555  0  0 -2  2 -2   -5.0    0.3
167.565  0  0 -2  2 -1    0.2    0.0
168.554  0 -1 -2  2 -2   -0.2    0.0
173.655  1  0  0 -2  0   -0.5    0.0
173.665  1  0  0 -2  1   -0.1    0.0
175.445 -1  0  0  0 -1    0.1    0.0
175.455 -1  0  0  0  0   -2.1    0.1
175.465 -1  0  0  0  1   -0.4    0.0
183.555  0  0  0 -2  0   -0.2    0.0
185.355 -2  0  0  0  0   -0.1    0.0
185.555  0  0 -2  0 -2   -0.6    0.0
185.565  0  0 -2  0 -1   -0.4    0.0
185.575  0  0 -2  0  0   -0.1    0.0
195.455 -1  0 -2  0 -2   -0.1    0.0
195.465 -1  0 -2  0 -1   -0.1    0.0
];
% dC21 dS21 computations
dC21 = 0;
dS21 = 0;
[sz1 sz2] = size(table65a);                                                 % [sz1 sz2] = size(table65b);
for i = 1 : sz1
    % thetaf (in radians)
    thetaf = 1 * (thetag + pi) - table65a(i,2:6) * [F1 F2 F3 F4 F5]';       % thetaf = 1 * (thetag + pi) - table65b(i,2:6) * [F1 F2 F3 F4 F5]';
    dC21 = dC21 + table65a(i,7) * 10^(-12) * sin(thetaf) + table65a(i,8) * 10^(-12) * cos(thetaf);
    dS21 = dS21 + table65a(i,7) * 10^(-12) * cos(thetaf) - table65a(i,8) * 10^(-12) * sin(thetaf);
    clear thetaf
end
clear sz1 sz2 i thetaf
dCnm(2+1,1+1) = dC21;
dSnm(2+1,1+1) = dS21;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dC22,dS22 corrections 
% Table 6.5c
% Matrix format : Columns :  Doodson_No. l l' F D Omega Amp.(ip)
table65c  = [
245.655  1  0  2  0  2   -0.3
255.555  0  0  2  0  2   -1.2
];
          
% dC22 dS22 computations
dC22 = 0;
dS22 = 0;
for i = 1 : 2
    % thetaf (in radians)
    thetaf = 2 * (thetag + pi) - table65c(i,2:6) * [F1 F2 F3 F4 F5]';
    dC22 = dC22 + table65c(i,7) * 10^(-12) * cos(thetaf); 
    dS22 = dS22 - table65c(i,7) * 10^(-12) * sin(thetaf); 
    clear thetaf
end
clear i thetaf
dCnm(2+1,2+1) = dC22;
dSnm(2+1,2+1) = dS22;
% dCnm(2+1,2+1) = 0;
% dSnm(2+1,2+1) = 0;

