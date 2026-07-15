%--------------------------------------------------------------------------
% GEORB models' data
%--------------------------------------------------------------------------

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
data_path = fullfile(pwd_path,'..','data');
cd(data_path);
data_path = pwd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s \n', 'GEORB scripts:');
fprintf('%s \n', 'Script for download input data :: Models data');
fprintf('%s%s \n\n', 'Download and Save data in ', data_path);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth Orientation Parameters by IERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s ', 'Earth Orientation Parameters :: EOP-C04 :: ');
url_read = 'https://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now';
save_fname = 'eopc04_IAU2000.62-now';
outpath = websave(save_fname , url_read);
fprintf('%s \n\n', 'downloaded');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leap Seconds 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s ', 'IERS Bulletin C Leap Seconds table:: IERS Leap_Second.dat :: ');
url_read = 'https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat';
save_fname = 'Leap_Second.dat';
outpath = websave(save_fname , url_read);
fprintf('%s \n\n', 'downloaded');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('%s%s \n\n', 'Data have been downloaded and saved in :: ', data_path);

% Direct to the central path of the package
georb_path = fullfile(pwd_path,'..');
cd(georb_path);
