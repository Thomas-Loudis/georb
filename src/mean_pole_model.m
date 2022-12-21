function [xp_mean, yp_mean] = mean_pole_model(mjd)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean_pole_model : Mean Pole model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Mean Pole model 
%
% Input arguments
% - mjd     : MJD including fraction of the day 
%
% Output arguments:
% - xp_mean    : Mean Pole coordiantes xp at input epoch in mas (milli-arcsec)
% - yp_mean    : Mean Pole coordiantes yp at input epoch in mas (milli-arcsec)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
%  Computed dCnm and dSnm are stored into lower triangular matrices.
%  Coefficient dCnm corresponds to matrix element dCnm(n+1,m+1).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             8 October 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients of the IERS (2010) mean pole model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 7.7
% Matrix format : Degree i,  xp(i) / mas yr^-i,  yp(i) / / mas yr^-i

% Coefficients until 2010
table77_meanpole_2010_pro  = [
0  55.974  346.346
1  1.8243  1.7896
2  0.18413  -0.10729
3  0.007024  -0.000908
];

% Coefficients after 2010
table77_meanpole_2010_post = [
0  23.513  358.891
1  7.6141  -0.6287
2  0  0
3  0  0
];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Year time argument
[sec,day,month,year]   = MJD_inv(mjd);
[JD_year_0,MJD_year_0] = MJD_date(0,1,1,year);
year_fr = (mjd - MJD_year_0) / 365.25;
t_year = year + year_fr;

if t_year > 2010.0
  meanpole_coeff_table = table77_meanpole_2010_post;  
else
  meanpole_coeff_table = table77_meanpole_2010_pro;      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean Pole coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0   = 2000.0;
xp_t = 0;
yp_t = 0;
[N_coeff, N_coll] = size(meanpole_coeff_table);
for k = 1 : N_coeff
    delta_t = t_year - t0;
    degree_i = meanpole_coeff_table(k,1);
    xp_i     = meanpole_coeff_table(k,2);
    yp_i     = meanpole_coeff_table(k,3);
    xp_t = xp_t + delta_t^degree_i * xp_i ;
    yp_t = yp_t + delta_t^degree_i * yp_i ;
end
xp_mean = xp_t;
yp_mean = yp_t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
