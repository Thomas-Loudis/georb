function [X_PN] = pn_series(table_series, arg_matrix) 

%-------------------------------------------------------------------------- 
% Copyright (C) 2007-present  Thomas (Loudis) Papanikolaou  
%  
% GEORB :: 'Gravity and Precise Orbit Determination system' 
%  
% GEORB is a software for Precise Orbit Determination (POD), gravity field  
% recovery and mission design 
%  
% GEORB is licensed under the GNU Affero General Public License v3.0 while  
% for information regarding dual-license options visit the official website 
% www.georb.gr  
%-------------------------------------------------------------------------- 
 

%-------------------------------------------------------------------------- 
% Function:  pn_series 
%-------------------------------------------------------------------------- 
% Purpose: 
%  Precession-Nutation model (IAU2006/2000A) non-polynomial series  
%  
%-------------------------------------------------------------------------- 
% Input arguments: 
% - MJD:       Epoch's Modified Julian Day in Terrestrial Time scale 
% 
% Output arguments: 
% - X_PN:      Coordinate series sum (X or Y) of the Precession-Nutation model in microarcsec 
% 
%-------------------------------------------------------------------------- 
% Thomas Loudis Papanikolaou                                   30 June 2026 
%-------------------------------------------------------------------------- 


%-------------------------------------------------------------------------- 
% PN model Series 
%-------------------------------------------------------------------------- 
% X series J=0
% table_series = X_series_j_0;
X_sum = 0;

% Fundamental arguments of Lunisolar and Planetary nutation arguments
arg_fund = arg_matrix;

%-------------------------------------------------------------------------- 
[n_series, m_col] = size(table_series);
for i_series = 1 : n_series
    % Fundamental Arguments :: LuniSolar and Planetary arguments 
    % arg_fund = [F1; F2; F3; F4; F5; F6; F7; F8; F9; F10; F11; F12; F13; F14];

    arguments_multipliers = table_series(i_series, 4 : 17);

    % Arguments product in radians 
    ARG_sum_i = arguments_multipliers * arg_fund; 

    % Xsum_J0 = Sum_i[a_{s,0})_i * sin(ARG) + a_{c,0})_i * cos(ARG)]
    a_s_j = table_series(i_series,2);
    a_c_j = table_series(i_series,3);

    X_i_j = a_s_j * sin(ARG_sum_i) + a_c_j * cos(ARG_sum_i);

    X_sum = X_sum + X_i_j;
end
%-------------------------------------------------------------------------- 

X_PN = X_sum;
%-------------------------------------------------------------------------- 

