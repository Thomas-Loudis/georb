function [GM,ae,Cnm,Snm,sCnm,sSnm,nmax, tide_system] = gfc_tv2(egmfilename,n_trunc,sigma_shc,mjd_t)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: gfc_tv2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read time-variable gravity model's spherical harmonic coefficients from
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
%       ae = sscanf(line_ith,'%*s %e %*')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static gravity coefficients
    test = strcmp(str1,'gfc');
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-variable gravity coefficients
    test = strcmp(str1,'gfct');	
    if test == 1
      n_ith = sscanf(line_ith,'%*s %d %*');      
      if (n_ith <= n_trunc) | (n_trunc == -1) 
          % Epoch test
          t0      = sscanf(line_ith,'%*s %*s %*s %*s %*s %*s %*s %s %*');          
          tn      = sscanf(line_ith,'%*s %*s %*s %*s %*s %*s %*s %*s %s %*');          
          % Epochs test
      	  Y_t0 = str2num(t0(1:4));  
      	  M_t0 = str2num(t0(5:6));  
      	  D_t0 = str2num(t0(7:8));
		  sec_t0 = 0;
		  [jd_t0,mjd_t0] = MJD_date(sec_t0,D_t0,M_t0,Y_t0);
		  
      	  Y_tn = str2num(tn(1:4));  
      	  M_tn = str2num(tn(5:6));  
      	  D_tn = str2num(tn(7:8));
		  sec_tn = 0;	  
		  [jd_tn,mjd_tn] = MJD_date(sec_tn,D_tn,M_tn,Y_tn);

		if (mjd_t > mjd_t0) && (mjd_t < mjd_tn)	                    
          % read SHC degree, order and respective coefficients values
          n = n_ith;
          m = sscanf(line_ith,'%*s %*d %d %*');

          % C,S at to              
          Cnm_t0  = str2num( sscanf(line_ith,'%*s %*s %*s %s %*') );
          Snm_t0  = str2num( sscanf(line_ith,'%*s %*s %*s %*s %s %*') );          
          sCnm_t0 = str2num( sscanf(line_ith,'%*s %*s %*s %*s %*s %s %*') );          
          sSnm_t0 = str2num( sscanf(line_ith,'%*s %*s %*s %*s %*s %*s %s %*') );          
          
          % Trend terms
          line_ith = fgetl(fid);
		  test_to = 0;
		  while test_to == 0              
          t0      = sscanf(line_ith,'%*s %*s %*s %*s %*s %*s %*s %s %*');          
          tn      = sscanf(line_ith,'%*s %*s %*s %*s %*s %*s %*s %*s %s %*');          
          % Epochs test
      	  Y_t0 = str2num(t0(1:4));  
      	  M_t0 = str2num(t0(5:6));  
      	  D_t0 = str2num(t0(7:8));
		  sec_t0 = 0;
		  [jd_t0,mjd_t0] = MJD_date(sec_t0,D_t0,M_t0,Y_t0);
		  
      	  Y_tn = str2num(tn(1:4));  
      	  M_tn = str2num(tn(5:6));  
      	  D_tn = str2num(tn(7:8));
		  sec_tn = 0;	  
		  [jd_tn,mjd_tn] = MJD_date(sec_tn,D_tn,M_tn,Y_tn);

		  if (mjd_t > mjd_t0) && (mjd_t < mjd_tn)
          test_to = 1;
		  % Trend
          Cnm_trend  = sscanf(line_ith,'%*s %*d %*d %f %*f %*f %*f %*');
          Snm_trend  = sscanf(line_ith,'%*s %*d %*d %*f %f %*f %*f %*');
          sCnm_trend = sscanf(line_ith,'%*s %*d %*d %*f %*f %f %*f %*');
          sSnm_trend = sscanf(line_ith,'%*s %*d %*d %*f %*f %*f %f %*');	  
          else
              line_ith = fgetl(fid);
          end
          end
          
          % annual and semi-annual terms
          line_ith = fgetl(fid);
		  test_to = 0;
		  while test_to == 0              
          t0      = sscanf(line_ith,'%*s %*s %*s %*s %*s %*s %*s %s %*');          
          tn      = sscanf(line_ith,'%*s %*s %*s %*s %*s %*s %*s %*s %s %*');          
          % Epochs test
      	  Y_t0 = str2num(t0(1:4));  
      	  M_t0 = str2num(t0(5:6));  
      	  D_t0 = str2num(t0(7:8));
		  sec_t0 = 0;
		  [jd_t0,mjd_t0] = MJD_date(sec_t0,D_t0,M_t0,Y_t0);
		  
      	  Y_tn = str2num(tn(1:4));  
      	  M_tn = str2num(tn(5:6));  
      	  D_tn = str2num(tn(7:8));
		  sec_tn = 0;	  
		  [jd_tn,mjd_tn] = MJD_date(sec_tn,D_tn,M_tn,Y_tn);

		  if (mjd_t > mjd_t0) && (mjd_t < mjd_tn)	
          test_to = 1;
            % annual terms		  
            % acos1
            Cnm_acos1  = sscanf(line_ith,'%*s %*d %*d %f %*f %*f %*f %*f %*f %*f %*');
            Snm_acos1  = sscanf(line_ith,'%*s %*d %*d %*f %f %*f %*f %*f %*f %*f %*');
            sCnm_acos1 = sscanf(line_ith,'%*s %*d %*d %*f %*f %f %*f %*f %*f %*f %*');
            sSnm_acos1 = sscanf(line_ith,'%*s %*d %*d %*f %*f %*f %f %*f %*f %*f %*');	  
            period1    = sscanf(line_ith,'%*s %*d %*d %*f %*f %*f %*f %*f %*f %f %*');	  				
            % asin1
            line_ith = fgetl(fid);
            Cnm_asin1  = sscanf(line_ith,'%*s %*d %*d %f %*f %*f %*f %*f %*f %*f %*');
            Snm_asin1  = sscanf(line_ith,'%*s %*d %*d %*f %f %*f %*f %*f %*f %*f %*');
            sCnm_asin1 = sscanf(line_ith,'%*s %*d %*d %*f %*f %f %*f %*f %*f %*f %*');
            sSnm_asin1 = sscanf(line_ith,'%*s %*d %*d %*f %*f %*f %f %*f %*f %*f %*');	  
            period1    = sscanf(line_ith,'%*s %*d %*d %*f %*f %*f %*f %*f %*f %f %*');	  

            % semi-annual terms		  
            % acos2
            line_ith = fgetl(fid);
            Cnm_acos2  = sscanf(line_ith,'%*s %*d %*d %f %*f %*f %*f %*f %*f %*f %*');
            Snm_acos2  = sscanf(line_ith,'%*s %*d %*d %*f %f %*f %*f %*f %*f %*f %*');
            sCnm_acos2 = sscanf(line_ith,'%*s %*d %*d %*f %*f %f %*f %*f %*f %*f %*');
            sSnm_acos2 = sscanf(line_ith,'%*s %*d %*d %*f %*f %*f %f %*f %*f %*f %*');	  
            period2    = sscanf(line_ith,'%*s %*d %*d %*f %*f %*f %*f %*f %*f %f %*');	  				
            % asin2
            line_ith = fgetl(fid);
            Cnm_asin2  = sscanf(line_ith,'%*s %*d %*d %f %*f %*f %*f %*f %*f %*f %*');
            Snm_asin2  = sscanf(line_ith,'%*s %*d %*d %*f %f %*f %*f %*f %*f %*f %*');
            sCnm_asin2 = sscanf(line_ith,'%*s %*d %*d %*f %*f %f %*f %*f %*f %*f %*');
            sSnm_asin2 = sscanf(line_ith,'%*s %*d %*d %*f %*f %*f %f %*f %*f %*f %*');	  
            period2    = sscanf(line_ith,'%*s %*d %*d %*f %*f %*f %*f %*f %*f %f %*');	              
          else
              line_ith = fgetl(fid);
          end
          end

		  test_to = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%G(t)=gfct(t0)+trnd*(t-t0)+asin1*sin(2pi/p1*(t-t0))+acos1*cos(2pi/p1*(t-t0))
%                         +asin2*sin(2pi/p2*(t-t0))+acos2*cos(2pi/p2*(t-t0))
T_days = period1 * 365.25;
delta_t = (mjd_t - mjd_t0) / T_days;			
Cnm_i = Cnm_t0 + Cnm_trend * delta_t + Cnm_acos1 * cos(2*pi/period1 * delta_t) + Cnm_asin1 * sin(2*pi/period1 * delta_t) ...
    								 + Cnm_acos2 * cos(2*pi/period2 * delta_t) + Cnm_asin2 * sin(2*pi/period2 * delta_t) ;                               
Snm_i = Snm_t0 + Snm_trend * delta_t + Snm_acos1 * cos(2*pi/period1 * delta_t) + Snm_asin1 * sin(2*pi/period1 * delta_t) ...
									 + Snm_acos2 * cos(2*pi/period2 * delta_t) + Snm_asin2 * sin(2*pi/period2 * delta_t) ;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          Cnm(n+1,m+1) = Cnm_i;
          Snm(n+1,m+1) = Snm_i;          
          if sigma_shc ~= 0
              sCnm(n+1,m+1) = sCnm_t0;
              sSnm(n+1,m+1) = sSnm_t0;
          end
          clear n m
                                
      end
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
