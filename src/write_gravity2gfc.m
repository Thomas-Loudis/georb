function [gravity_model_filename] = write_gravity2gfc(gravity_model_filename, GM,radius,Cnm,Snm, sigma_Cnm, sigma_Snm,n_trunc,tide_system, time_period_of_data)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  write_gravity2gfc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Write gravity field solution into gravity models' gfc format 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - out_filename : Gravtiy model output file name (*.gfc)
% - GM           : Earth gravity constant  (m^3/sec^2)
% - ae           : Radius  (meters)
% - Cnm, Snm     : Spherical harmonic coefficients matrices
% - sCnm, sSnm   : Covariance matrix of sherical harmonic coefficients
% - n_trunc      : Degree limit for writing .gfc file
%                  (set to -1 for writing all spherical harmonic coefficients)
% - tide_system  : Tides system of the gravity model: 'tide_free' or 'zero_tide'
%
% Output arguments:
% - gravity_model_filename : Gravity field model output file name 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                               9 June  2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n_char = length(gravity_model_filename);
gravity_model_name = gravity_model_filename(1 : n_char-4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity model's Maximum Degree max_degree
[n, m] = size(Cnm);
Nmax = n - 1;
if n_trunc == -1
    max_degree = Nmax;
else
    max_degree = n_trunc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open file for writing
fid = fopen(gravity_model_filename,'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Header :: Gravity model general comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text_line = 'GEORB Gravity Field solution :: Unconstrained GRACE Follow-On gravity field solution';
fprintf(fid,'%s\n',text_line);
%fprintf(fid,'%s\n','');
fprintf(fid,'%s\n',' ');

text_line = 'Reference for citing when using the data:';
fprintf(fid,'%s\n',text_line);
fprintf(fid,'%s\n',' ');
text_line = 'Thomas Loudis Papanikolaou (2023). GEORB: Release for precise orbit determination of low Earth orbiters and satellite gravity missions, ';
fprintf(fid,'%s\n',text_line);
text_line = 'Software Impacts, doi: https://doi.org/10.1016/j.simpa.2023.100502';
fprintf(fid,'%s\n',text_line);

fprintf(fid,'%s\n',' ');

param_keyword = 'time_period_of_data';
fprintf(fid,'%s %s ',param_keyword, ' ');
fprintf(fid,'%s\n','');

fprintf(fid,'%s\n',' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Header :: Gravity model parameters of icegem format 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text_line = 'begin_of_head =================================================================================';
fprintf(fid,'%s\n',text_line);

param_keyword = 'modelname';
fprintf(fid,'%-23s %s ',param_keyword, gravity_model_name);
fprintf(fid,'%s\n','');

param_keyword = 'product_type';
fprintf(fid,'%-23s %s ',param_keyword, 'gravity_field');
fprintf(fid,'%s\n','');

param_keyword = 'earth_gravity_constant';
fprintf(fid,'%-23s %.10e ',param_keyword, GM);
fprintf(fid,'%s\n','');

param_keyword = 'radius';
fprintf(fid,'%-23s %.10e ',param_keyword, radius);
fprintf(fid,'%s\n','');

param_keyword = 'max_degree';
fprintf(fid,'%-23s %d ',param_keyword, max_degree);
fprintf(fid,'%s\n','');

param_keyword = 'norm';
fprintf(fid,'%-23s %s ',param_keyword, 'fully_normalized');
fprintf(fid,'%s\n','');

param_keyword = 'tide_system';
fprintf(fid,'%-23s %s ',param_keyword, tide_system);
fprintf(fid,'%s\n','');

param_keyword = 'errors';
fprintf(fid,'%-23s %s ',param_keyword, 'formal');
fprintf(fid,'%s\n','');

text_line = 'key      L    M         C                   S                sigma C             sigma S';
fprintf(fid,'%s\n',text_line);

text_line = 'end_of_head ===================================================================================';
fprintf(fid,'%s\n',text_line);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write sherical harmonic coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 0 : max_degree
    for m = 0 : n
c_coef = Cnm(n+1,m+1);
s_coef = Snm(n+1,m+1);
sigma_c = sigma_Cnm(n+1,m+1);
sigma_s = sigma_Snm(n+1,m+1);
fprintf(fid,'%-5s%5d%5d%20.12e%20.12e%20.12e%20.12e','gfc', n, m, c_coef, s_coef, sigma_c, sigma_s);
fprintf(fid,'%s\n',' ');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Close file
fclose(fid);

