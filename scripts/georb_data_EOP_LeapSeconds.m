%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB models' data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% georb_data_EOP_LeapSeconds.m script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  georb_data_EOP_LeapSeconds is a script file for downloading/updating 
%  the Earth Orientation Parameters (EOP) and Leap Seconds data 
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
data_output_isfolder = isfolder(data_path);
if data_output_isfolder == 0
[status, message, messageid] = mkdir(data_path);
end
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
