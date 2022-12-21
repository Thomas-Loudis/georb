%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB models' data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% georb_data_models.m script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  georb_data_models is a script file for downloading the basic models'
%  data required for the operation of the GEORB software  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                   18 June 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
format long e
fclose('all');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Folders path for input data (models)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pwd_path = pwd;
data_path_fname = '../data/';
data_path = fullfile(pwd_path,data_path_fname);
cd(data_path);
data_path = pwd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s \n', 'GEORB scripts:');
fprintf('%s \n', 'Script for download input data :: Models data');
fprintf('%s%s \n\n', 'Download and Save data in ', data_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field models by ICGEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GOCO06s
fprintf('%s', 'Gravity Field model :: GOCO06s :: ');
url_read   = 'http://icgem.gfz-potsdam.de/getmodel/gfc/32ec2884630a02670476f752d2a2bf1c395d8c8d6d768090ed95b4f04b0d5863/GOCO06s.gfc';
save_fname = 'GOCO06s.gfc';
outpath = websave(save_fname , url_read);
fprintf('%s \n', 'downloaded');

% EIGEN-6S4
fprintf('%s', 'Gravity Field model :: EIGEN-6S4 :: ');
url_read   = 'http://icgem.gfz-potsdam.de/getmodel/gfc/4cc119d62d3f83adce857914bcabdfa7ca2f91677ed2905934cf52698584b4c9/EIGEN-6S4%20(v2).gfc';
save_fname = 'EIGEN-6S4.gfc';
outpath = websave(save_fname , url_read);
fprintf('%s \n', 'downloaded');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ocean Tides model by CNES/GRGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s', 'Ocean Tides model :: FES2014b-v1 :: ');
url_read   = 'http://gravitegrace.get.obs-mip.fr/geofluid/fes2014b.v1.Cnm-Snm_Om1+Om2C20_with_S1_wave.POT_format.txt';
save_fname = 'fes2014b_v1.txt';
outpath = websave(save_fname , url_read);
fprintf('%s \n', 'downloaded');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Planetary and Lunar ephemeris by JPL/NASA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s', 'Planetary/Lunar ephemeris :: DE423 :: ');
url_read   = 'https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de423/ascp2000.423';
save_fname = 'ascp2000.423';
outpath = websave(save_fname , url_read);

url_read   = 'https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de423/header.423';
save_fname = 'header.423';
outpath = websave(save_fname , url_read);

fprintf('%s \n', 'downloaded');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Orientation Parameters by IERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s ', 'Earth Orientation Parameters :: EOP-C04 :: ');
url_read = 'https://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now';
save_fname = 'eopc04_IAU2000.62-now';
outpath = websave(save_fname , url_read);
fprintf('%s \n', 'downloaded');

fprintf('%s ', 'IERS Bulletin C Leap Seconds table:: IERS Leap_Second.dat :: ');
url_read = 'https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat';
save_fname = 'Leap_Second.dat';
outpath = websave(save_fname , url_read);
fprintf('%s \n', 'downloaded');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s%s \n\n', 'Data have been downloaded and saved in :: ', data_path);

% Direct to the central path of the package
gsynth_path = fullfile(pwd_path,'/../');
cd(gsynth_path);
