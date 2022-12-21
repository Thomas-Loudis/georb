function [otides_struct, delaunay_doodson_multipliers, dCnm_plus, dSnm_plus, dCnm_minus, dSnm_minus] = read_oceantides(fesmodel_filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_oceantides : Read Ocean Tides Model of FES series format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read data of ocean tide models from FES model series (FES2004, FES2014b) 
%
% Input arguments
% - fesmodel_filename : FES model e.g. FES2004, FES2014b, data file name 
%                       with the geopotential Cnm, Snm coefficients 
%
% Output arguments:
% - DelaunayNf_fes2004 : Delaunay variables multipliers for FES2004 tidal
%                        waves 
% - dCnm_plus  : FES2004 coefficients to be used in Eq.6... IERS Conv.2010
%   dSnm_plus    Coefficients are stored in these four 3D arrays
%   dCnm_minus
%   dSnm_minus
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Papanikolaou                                            June  2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
%  23/06/2012  Speed up reading of FES2004 model's coefficients and store
%              them in 3D arrays
%  28/06/2012  Revised according to the IERS Conventions 2010 updates at 
%              23/09/2011 & 14/10/2011
%              fes2004_Cnm-Snm.dat : revised file
%              Coefficients Unit is now 10^-11 and Coefficients corrections
% 
% 01/02/2019  Dr. Thomas Papanikolaou
%             Delaunay variables multipliers update: M4 tidal wave added 
% 01/12/2022  Dr. Thomas Loudis Papanikolaou
%             Major upgrade: Doodson argument number through calling the 
%             new funtion doodson_number.m; 
%             Delaunay variables multipliers has been replaced to dynamic
%             array;
%             Doodson and Delaunay multipliers are computed for each
%             model's tidal frequency and saved to structure array;
%             Function rename from tides_fes2004.m to read_oceantides.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients Unit : 10^-12 (IERS Conventions 2010)
% Coefficients Unit : 10^-11 (IERS Conventions 2010 updated 23/09/11)
fes2004_cfunit = 10^(-11);
%fprintf('%s %.0e \n','FES2004 Coefficients Unit :',fes2004_cfunit);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read FES model data file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(fesmodel_filename);
% linei = 0;
headeri = 1;
while headeri == 1
    HDlineread = fgetl(fid);
    str_1 = sscanf(HDlineread,'%s%*');
    test = strcmp(str_1,'Doodson');
    if test == 1
        headeri = 0;
    end
end

doodson_number_current = '0';
i_doodson_number = 0;
Nmax  = 0;
Nread = 0;
while (~feof(fid))
    lineread = fgetl(fid);
    str_1 = sscanf(lineread,'%s%*');    
    % Doodson number    
    doodson_number_char = str_1;
    test = strcmp(doodson_number_char, doodson_number_current);
    if test == 0
        i_doodson_number = i_doodson_number + 1;
        doodson_number_current = doodson_number_char;
        if Nread > Nmax
            Nmax = Nread;
        end
    end
    % Read and store d/o coefficients
    dCnmSnm = sscanf(lineread,'%*f %*s %d %d %f %f %f %f');
    Nread = dCnmSnm(1,1);
    Mread = dCnmSnm(2,1);
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation
Nmax_otides  = Nmax;
%Nmax_otides  = Nread
Nfreq_otides = i_doodson_number;
delaunay_doodson_multipliers = zeros(Nfreq_otides, 12);
dCnm_plus  = zeros(Nmax_otides+1,Nmax_otides+1,Nfreq_otides);
dSnm_plus  = zeros(Nmax_otides+1,Nmax_otides+1,Nfreq_otides);
dCnm_minus = zeros(Nmax_otides+1,Nmax_otides+1,Nfreq_otides);
dSnm_minus = zeros(Nmax_otides+1,Nmax_otides+1,Nfreq_otides);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read FES model data file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(fesmodel_filename);
% linei = 0;
headeri = 1;
while headeri == 1
    HDlineread = fgetl(fid);
%     linei = linei + 1;
    str_1 = sscanf(HDlineread,'%s%*');
    test = strcmp(str_1,'Doodson');
    if test == 1
        headeri = 0;
    end
end

doodson_number_current = '0';
i_doodson_number = 0;
while (~feof(fid))
    lineread = fgetl(fid);
%     linei = linei + 1;
    str_1 = sscanf(lineread,'%s%*');
    
    % Doodson number    
    doodson_number_char = str_1;
    test = strcmp(doodson_number_char,doodson_number_current);
    if test == 0
        i_doodson_number = i_doodson_number + 1;
        doodson_number_current = doodson_number_char;

        doodson_number_num = sscanf(lineread,'%f%*');
        % Delaunay and Doodson multipliers based on Doodson number
        [doodson_multipliers, delaunay_multipliers] = doodson_number(doodson_number_char);
        % Array of multipliers
        delaunay_doodson_multipliers (i_doodson_number,:) = [doodson_number_num delaunay_multipliers doodson_multipliers];
    end
    
    % Read and store d/o coefficients
    dCnmSnm = sscanf(lineread,'%*f %*s %d %d %f %f %f %f');
    Nread = dCnmSnm(1,1);
    Mread = dCnmSnm(2,1);
        
    % 3D arrays  ixjxk : n x m x frq
    dCnm_plus (Nread+1, Mread+1, i_doodson_number) = fes2004_cfunit * dCnmSnm(3,1);
    dSnm_plus (Nread+1, Mread+1, i_doodson_number) = fes2004_cfunit * dCnmSnm(4,1);
    dCnm_minus(Nread+1, Mread+1, i_doodson_number) = fes2004_cfunit * dCnmSnm(5,1);
    dSnm_minus(Nread+1, Mread+1, i_doodson_number) = fes2004_cfunit * dCnmSnm(6,1);
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strucutre array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i_struct = 0;
% 
% i_struct = i_struct + 1;
% otides_struct(i_struct).names  = 'delaunay_doodson_multipliers';
% otides_struct(i_struct).values = delaunay_doodson_multipliers;
% 
% i_struct = i_struct + 1;
% otides_struct(i_struct).names  = 'dCnm_plus';
% otides_struct(i_struct).values = dCnm_plus;
% 
% i_struct = i_struct + 1;
% otides_struct(i_struct).names  = 'dSnm_plus';
% otides_struct(i_struct).values = dSnm_plus;
% 
% i_struct = i_struct + 1;
% otides_struct(i_struct).names  = 'dCnm_minus';
% otides_struct(i_struct).values = dCnm_minus;
% 
% i_struct = i_struct + 1;
% otides_struct(i_struct).names  = 'dSnm_minus';
% otides_struct(i_struct).values = dSnm_minus;

otides_struct.degree = Nmax_otides;
otides_struct.delaunay_doodson_multipl = delaunay_doodson_multipliers;
otides_struct.dCnm_plus  = dCnm_plus;
otides_struct.dSnm_plus  = dSnm_plus;
otides_struct.dCnm_minus = dCnm_minus;
otides_struct.dSnm_minus = dSnm_minus;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
