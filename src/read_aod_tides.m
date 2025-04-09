function [GM,ae,nmax, dCnm_plus, dSnm_plus, dCnm_minus, dSnm_minus] = read_aod_tides(aod_tides_filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function :: read_aod_tides 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read the (atmospheric) tides data from the Atmosphere and Ocean
%  De-Aliasing (AOD) Level-1B (AOD1B) product
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - aod_tides_filename: AOD1B Tides filename
% 
% Output arguments:
% - GM:                 Earth gravity constant  (m^3/sec^2)
% - ae:                 radius  (meters)
% - nmax:               Cnm and Snm matrices maximum degree
% - Cnm, Snm:           normalized spherical harmonics coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
% Coefficient Cnm corresponds to matrix element Cnm(n+1,m+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                            2 December 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(aod_tides_filename);
% while (~feof(fid))
header_status = 1;
while (header_status == 1)
    line_ith = fgetl(fid);
    str_test = sscanf(line_ith,'%s%1c%s %*');

    test = strcmp(str_test,'MAXIMUM DEGREE');
    if test == 1
      Nmax_aod = sscanf(line_ith,'%*s %*s %*s %d %*');
    end
    
    test = strcmp(str_test,'CONSTANT GM');
    if test == 1
      GM = sscanf(line_ith,'%*s %*s %*s %*s %f %*');
    end

    test = strcmp(str_test,'CONSTANT A');
    if test == 1
      constant_a = sscanf(line_ith,'%*s %*s %*s %*s %f %*');
    end
    
    test = strcmp(str_test,'CONSTANT FLAT');
    if test == 1
      constant_flat = sscanf(line_ith,'%*s %*s %*s %*s %f %*');
    end

    test = strcmp(str_test,'CONSTANT OMEGA');
    if test == 1
      constant_omega = sscanf(line_ith,'%*s %*s %*s %*s %f %*');
    end

    keyword_test = sscanf(line_ith,'%s%1c%s%1c%s%1c%s %*');
    test = strcmp(keyword_test,'NUMBER OF DATA SETS');
    if test == 1
      data_sets_number = sscanf(line_ith,'%*s %*s %*s %*s %*s %d %*');
    end
    
    keyword_test = sscanf(line_ith,'%s%1c%s%1c%s%1c%s %*');
    test = strcmp(keyword_test,'NUMBER OF DATA RECORDS');
    if test == 1
      data_rec_number = sscanf(line_ith,'%*s %*s %*s %*s %*s %d %*');
    end  
    
    str_test = sscanf(line_ith,'%s%1c%s%1c%s %*');
    test = strcmp(str_test,'END OF HEADER');
    if test == 1
      header_status = 0;
    end        
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ae   = constant_a;
nmax = Nmax_aod;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arrays initialization
dCnm_plus  = zeros(Nmax_aod+1,Nmax_aod+1);
dCnm_minus = zeros(Nmax_aod+1,Nmax_aod+1);

dSnm_plus  = zeros(Nmax_aod+1,Nmax_aod+1);
dSnm_minus = zeros(Nmax_aod+1,Nmax_aod+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_rec_number_cos = data_rec_number / 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Stokes coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(aod_tides_filename);
header_status = 1;
while (header_status == 1)
    line_ith = fgetl(fid);
    str_test = sscanf(line_ith,'%s%1c%s%1c%s %*');
    test = strcmp(str_test,'END OF HEADER');
    if test == 1
      header_status = 0;
    end    
end

% DATA SET 01:  16471 COEFFICIENTS OF TYPE cos
line_ith = fgetl(fid);
keyword_data_set = sscanf(line_ith,'%s%1c%s%1c%s %*');
N_coefficients = sscanf(line_ith,'%*s%*s%*s %d %*');

i_datarec_cos = 1;
while (i_datarec_cos <= N_coefficients)
    line_ith = fgetl(fid);
    keyword_data_set = sscanf(line_ith,'%s%1c%s %*');    
    % Read and store Stokes coeffiecients of all data sets effects
    % n_i   = sscanf(line_ith,'%d %*')
    % m_i   = sscanf(line_ith,'%*d %d %*')
    % Cnm_plus  = sscanf(line_ith,'%*d %*d %f %*');
    % Cnm_minus = sscanf(line_ith,'%*d %*d %*f %f %*');
    [data_ith_vec] = sscanf(line_ith,'%d %d %f %f %*');
    n_i = data_ith_vec(1,1);
    m_i = data_ith_vec(2,1);
    Cnm_plus = data_ith_vec(3,1);
    Cnm_minus = data_ith_vec(4,1);      
    % Store coefficients to overall 3dimensional matrices
    dCnm_plus (n_i + 1, m_i + 1) = Cnm_plus;
    dCnm_minus(n_i + 1, m_i + 1) = Cnm_minus;
    % Next data record
    i_datarec_cos = i_datarec_cos + 1;
end

% DATA SET 02:  16471 COEFFICIENTS OF TYPE sin
line_ith = fgetl(fid);
keyword_data_set = sscanf(line_ith,'%s%1c%s%1c%s %*');
N_coefficients = sscanf(line_ith,'%*s%*s%*s %d %*');

i_datarec_cos = 1;
while (i_datarec_cos <= N_coefficients)
    line_ith = fgetl(fid);
    keyword_data_set = sscanf(line_ith,'%s%1c%s %*');    
    % Read and store Stokes coeffiecients of all data sets effects
    % n_i   = sscanf(line_ith,'%d %*');
    % m_i   = sscanf(line_ith,'%*d %d %*');
    % Snm_plus  = sscanf(line_ith,'%*d %*d %f %*');
    % Snm_minus = sscanf(line_ith,'%*d %*d %*f %f %*');    
    [data_ith_vec] = sscanf(line_ith,'%d %d %f %f %*');
    n_i = data_ith_vec(1,1);
    m_i = data_ith_vec(2,1);
    Snm_plus = data_ith_vec(3,1);
    Snm_minus = data_ith_vec(4,1);          
    % Store coefficients to overall 3dimensional matrices
    dSnm_plus (n_i + 1, m_i + 1) = Snm_plus;
    dSnm_minus(n_i + 1, m_i + 1) = Snm_minus;
    % Next data record
    i_datarec_cos = i_datarec_cos + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
