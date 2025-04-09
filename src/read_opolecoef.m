function [Anm_R,Anm_I,Bnm_R,Bnm_I,n_max] = read_opolecoef(coeffilename,n_trunc)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_opolecoef: Read Desai Ocean Pole tide coefficients Anm, Bnm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read Desai Ocean Pole tide coefficients Anm, Bnm from the data file
%  provided by IERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - coeffilename:       Coefficients data file name
% - n_trunc:            Degree truncation limit for reading coefficients 
%                       (set to -1 for reading all the SHC without limit)
% Output arguments:
% - Anm_R:              A coefficients Real part
% - Anm_I:              A coefficients Imaginary part
% - Bnm_R:              B coefficients Real part
% - Bnm_I:              B coefficients Imaginary part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
% Coefficient Cnm corresponds to matrix element Cnm(n+1,m+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou, DORUS Space Lab            31 August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Desai coefficients from IERS file  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(coeffilename);
iline = 0;
while (~feof(fid))
    line_ith = fgetl(fid);
    iline = iline + 1;
    if iline > 1
      n_ith = sscanf(line_ith,'%d %*');
    end
end
fclose(fid);
n_max = n_ith;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arrays initialization
if n_trunc == -1
    n_sz = n_max;
    m_sz = n_max;
else    
    n_sz = n_trunc;
    m_sz = n_trunc;
end
Anm_R = zeros(n_sz+1,n_sz+1);
Anm_I = zeros(n_sz+1,n_sz+1);
Bnm_R = zeros(n_sz+1,n_sz+1);
Bnm_I = zeros(n_sz+1,n_sz+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(coeffilename);
iline = 0;
while (~feof(fid))
    line_ith = fgetl(fid);
    iline = iline + 1;
    if iline > 1
      data_ith_vec = sscanf(line_ith,'%d %d %e %e %e %e %*');
      n              = data_ith_vec(1,1);
      m              = data_ith_vec(2,1);
      Anm_R(n+1,m+1) = data_ith_vec(3,1);
      Bnm_R(n+1,m+1) = data_ith_vec(4,1);
      Anm_I(n+1,m+1) = data_ith_vec(5,1);
      Bnm_I(n+1,m+1) = data_ith_vec(6,1);      
    end
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
