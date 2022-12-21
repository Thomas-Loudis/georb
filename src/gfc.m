function [GM,ae,Cnm,Snm,sCnm,sSnm,nmax,tide_system] = gfc(egmfilename,n_trunc,sigma_shc)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gfc: Gravity Field Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read gravity model's spherical harmonic coefficients from
%  the .gfc file in gfc format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - egmfilename:        Gravtiy model .gfc file's name
% - n_trunc:            Degree limit for reading .gfc file
%                       (set to -1 for reading all the SHC without limit)
% - sigma_shc:          Set to any value or to 0 value for reading or not
%                       the sigmaC and sigmaS from .gfc file
% Output arguments:
% - GM:                 Earth gravity constant  (m^3/sec^2)
% - ae:                 radius  (meters)
% - Cnm, Snm:           normalized spherical harmonics coefficients
% - sCnm,sSnm:          errors, coefficients covariances
% - nmax:               Cnm and Snm matrices maximum degree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
% Coefficient Cnm corresponds to matrix element Cnm(n+1,m+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                    June 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 31/03/2013   Revised for minimizing CPU time for higher degrees (e.g.
%              1000, 2000) by initializing Cnm,Snm matrices
% 25/08/2013   Upgrade for reading efficiently the ICGEM format and solve
%              gfc reading problems for the case of GOCE TIM models format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read .gfc file and format Cnm,Snm and sCnm,sSnm lower triangular matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(egmfilename);
while (~feof(fid))
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');

    test = strcmp(str1,'earth_gravity_constant');
    if test == 1
%       GM = sscanf(line_ith,'%*s %e %*');
      GM = str2num( sscanf(line_ith,'%*s %s %*') );
    end
    clear test

    test = strcmp(str1,'radius');
    if test == 1
%       ae = sscanf(line_ith,'%*s %f %*')
      ae = str2num( sscanf(line_ith,'%*s %s %*') );
    end
    clear test

    test = strcmp(str1,'max_degree');
    if test == 1
      Nmax_egm = sscanf(line_ith,'%*s %d %*');
    end
    clear test

    test = strcmp(str1,'tide_system');
    if test == 1
      tide_system = sscanf(line_ith,'%*s %s %*');
    end
    clear test    
end
fclose(fid);
clear fid line_ith str1 test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arrays initialization
if n_trunc == -1
    n_sz = Nmax_egm;
    m_sz = Nmax_egm;
else    
    n_sz = n_trunc;
    m_sz = n_trunc;
end
Cnm = zeros(n_sz+1,n_sz+1);
Snm = zeros(n_sz+1,n_sz+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(egmfilename);
while (~feof(fid))
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');

    test = strcmp(str1,'gfc');
    if test == 0
        test = strcmp(str1,'gfct');
    end
    if test == 1
      n_ith = sscanf(line_ith,'%*s %d %*');
      if (n_ith <= n_trunc) | (n_trunc == -1)
          % read SHC degree, order and respective coefficients values
          n = n_ith;
          m = sscanf(line_ith,'%*s %*d %d %*');
%           Cnm(n+1,m+1) = sscanf(line_ith,'%*s %*d %*d %f %*');
%           Snm(n+1,m+1) = sscanf(line_ith,'%*s %*d %*d %*f %f %*');
          Cnm(n+1,m+1) = str2num( sscanf(line_ith,'%*s %*s %*s %s %*') );
          Snm(n+1,m+1) = str2num( sscanf(line_ith,'%*s %*s %*s %*s %s %*') );          
          if sigma_shc ~= 0
%               sCnm(n+1,m+1) = sscanf(line_ith,'%*s %*d %*d %*f %*f %f %*');
%               sSnm(n+1,m+1) = sscanf(line_ith,'%*s %*d %*d %*f %*f %*f %f %*');
              sCnm(n+1,m+1) = str2num( sscanf(line_ith,'%*s %*s %*s %*s %*s %s %*') );
              sSnm(n+1,m+1) = str2num( sscanf(line_ith,'%*s %*s %*s %*s %*s %*s %s %*') );
          end
          clear n m
      end
    end
    clear test str1  
end
fclose(fid);
clear fid line_ith str1 test 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write zero to values of the coefficients C10,C11,S10,S11 and the
% corresponding sigmaC,sigmaS.
[sz1,sz2]=size(Cnm);
if n_trunc == 1 && sz1 == 1
    Cnm(1+1,1+1) = 0;
    Snm(1+1,1+1) = 0;
    if sigma_shc ~= 0
        sCnm(1+1,1+1) = 0;
        sSnm(1+1,1+1) = 0;
    end
end
clear sz1 sz2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sigma_shc == 0
    sCnm = 0;
    sSnm = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maximum degree (n) and order (m)
[nmax n2] = size(Cnm);
clear n2
nmax = nmax-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
