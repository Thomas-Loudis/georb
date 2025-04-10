function [GM,ae,Cnm,Snm,sCnm,sSnm,nmax, tide_system] = gfc_tv1(egmfilename,n_trunc,sigma_shc,mjd_t)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: gfc_tv1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read time-vraiable gravity model's spherical harmonic coefficients from
%  the .gfc file in gfc format.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - egmfilename:        EGM .gfc file's name
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
% while (~feof(fid))
header_status = 1;
while (header_status == 1)  
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');

    test = strcmp(str1,'earth_gravity_constant');
    if test == 1
      GM = sscanf(line_ith,'%*s %e %*');
    end

    test = strcmp(str1,'radius');
    if test == 1
      ae = sscanf(line_ith,'%*s %e %*');
    end

    test = strcmp(str1,'max_degree');
    if test == 1
      Nmax_egm = sscanf(line_ith,'%*s %d %*');
    end

    test = strcmp(str1,'tide_system');
    if test == 1
      tide_system = sscanf(line_ith,'%*s %s %*');
    end

    test = strcmp(str1,'end_of_head');
    if test == 1
      header_status = 0;
    end    

end
fclose(fid);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static gravity coefficients
    test = strcmp(str1,'gfc');
    if test == 1
      data_ith_vec = sscanf(line_ith,'%*s %d %d %e %e %e %e %*');
      n_ith = data_ith_vec(1,1);
      if (n_ith <= n_trunc) | (n_trunc == -1)
          % read SHC degree, order and respective coefficients values
          n = n_ith;
          m = data_ith_vec(2,1);
          Cnm(n+1,m+1) = data_ith_vec(3,1);
          Snm(n+1,m+1) = data_ith_vec(4,1);         
          if sigma_shc ~= 0
              sCnm(n+1,m+1) = data_ith_vec(5,1);
              sSnm(n+1,m+1) = data_ith_vec(6,1);         
          end
      end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-variable gravity coefficients
    test = strcmp(str1,'gfct');	
    if test == 1        
      data_ith_vec = sscanf(line_ith,'%*s %d %d %e %e %e %e %e %*');
      % n_ith = sscanf(line_ith,'%*s %d %*');
      n_ith = data_ith_vec(1,1);
      if (n_ith <= n_trunc) | (n_trunc == -1)
          % read SHC degree, order and respective coefficients values
          n = n_ith;
          t0      = sscanf(line_ith,'%*s %*s %*s %*s %*s %*s %*s %s %*');          		  
          m       = data_ith_vec(2,1);
          Cnm_t0  = data_ith_vec(3,1);
          Snm_t0  = data_ith_vec(4,1);          
          sCnm_t0 = data_ith_vec(5,1);          
          sSnm_t0 = data_ith_vec(6,1);          

          % Read next (3) lines for "trend,acos,asin" coefficients  
		  % Trend
		  line_ith = fgetl(fid);
          data_ith_vec = sscanf(line_ith,'%*s %d %d %e %e %e %e %e %*');
          Cnm_trend  = data_ith_vec(3,1);
          Snm_trend  = data_ith_vec(4,1);          
          sCnm_trend = data_ith_vec(5,1);          
          sSnm_trend = data_ith_vec(6,1);          
         
		  % acos
		  line_ith = fgetl(fid);
          data_ith_vec = sscanf(line_ith,'%*s %d %d %e %e %e %e %e %*');
          Cnm_acos  = data_ith_vec(3,1);
          Snm_acos  = data_ith_vec(4,1);          
          sCnm_acos = data_ith_vec(5,1);          
          sSnm_acos = data_ith_vec(6,1);          
          period    = data_ith_vec(7,1);          

		  % asin
		  line_ith = fgetl(fid);
          data_ith_vec = sscanf(line_ith,'%*s %d %d %e %e %e %e %e %*');
          Cnm_asin  = data_ith_vec(3,1);
          Snm_asin  = data_ith_vec(4,1);          
          sCnm_asin = data_ith_vec(5,1);          
          sSnm_asin = data_ith_vec(6,1);          
          period    = data_ith_vec(7,1);            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	  
		  % delta_t = (t_ys - t0_ys) / period
		  T_days = period * 365.25;
      	  Y_t0 = sscanf(t0(1:4),'%d%*');  
      	  M_t0 = sscanf(t0(5:6),'%d%*'); 
      	  D_t0 = sscanf(t0(7:8),'%d%*');
		  sec_t0 = 0;	  
		  [jd_t0,mjd_t0] = MJD_date(sec_t0,D_t0,M_t0,Y_t0);
		  delta_t = (mjd_t - mjd_t0) / T_days;			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		  % V(t) = gfct + trnd * (t-t0)/T + acos*cos( 2pi * (t-t0)/T ) + asin*sin( 2pi * (t-t0)/T )
Cnm_i = Cnm_t0 + Cnm_trend * delta_t + Cnm_acos * cos(2*pi * delta_t) + Cnm_asin * sin(2*pi * delta_t);
Snm_i = Snm_t0 + Snm_trend * delta_t + Snm_acos * cos(2*pi * delta_t) + Snm_asin * sin(2*pi * delta_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          Cnm(n+1,m+1) = Cnm_i;
          Snm(n+1,m+1) = Snm_i;          
          if sigma_shc ~= 0
              sCnm(n+1,m+1) = sCnm_t0;
              sSnm(n+1,m+1) = sSnm_t0;
          end
      end
    end
end
fclose(fid);
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
nmax = nmax-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
