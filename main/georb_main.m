%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 GEORB                                   %
%                                                                         %
%                                                                         %
%              Gravity and precisE ORBit determination system             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB is a software for Precise Orbit Determination (POD) of Low Earth
% Orbiters (LEOs), gravity Field modelling, data analysis of satellite
% gravity missions and design of future space missions. 
% 
% GEORB current version supports the orbit determination and accelerometer
% calibration modelling of the missions:
% - Gravity Recovery And Climate Experiment (GRACE) mission 
% - GRACE Follow-On mission  
%
% The source code has been written by Thomas since 2007. 
% GEORB was released as open source in 2022 and is available through the
% repository https://github.com/Thomas-Loudis/georb 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                                  July 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Citing GEORB:
%
% Papanikolaou T. (2022). Precise orbit determination and accelerometer
% data modelling of the GRACE Follow-On mission, GRACE/GRACE-FO Science
% Team Meeting 2022, Potsdam, Germany, 18–20 Oct 2022, GSTM2022-90,
% https://doi.org/10.5194/gstm2022-90, 2022.   
%
% References:
%
% Papanikolaou T. (2012). Dynamic modelling of satellite orbits in the 
% frame of contemporary space geodesy missions, Ph.D. Dissertation, 
% Aristotle University of Thessaloniki (AUTH), Greece.
%
% Papanikolaou T., Tsoulis D. (2016). Assessment of numerical integration
% methods in the context of low Earth orbits and inter-satellite
% observation analysis, Acta Geodetica et Geophysica, Vol. 51, No. 4, pp.
% 619-641, doi: 10.1007/s40328-016-0159-3.   
%
% Papanikolaou T., Tsoulis D. (2018). Assessment of Earth gravity field
% models in the medium to high frequency spectrum based on GRACE and GOCE
% dynamic orbit analysis, Geosciences, 8(12):441, doi:
% 10.3390/geosciences8120441 
%
% Papanikolaou T. (2022). Precise Orbit Determination and accelerometry
% calibration modelling of the GRACE Follow-On mission, Nordic Geodetic
% Commission (NKG) General Assembly, 5-8 September 2022, Copenhagen,
% Denmark.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% georb_main.m script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  georb_main is the main script file that calls the GEORB source code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 05/07/2022, Thomas Loudis Papanikolaou 
%             General upgrade of the source code for supporting the release
%             as open source software
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
format long e
fclose('all');
% Restore to default path 
restoredefaultpath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Folders path for source code, configuration files, input data, output
% data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
main_path = pwd;
%addpath(genpath(main_path));

data_path = fullfile(main_path,'/../data/');
addpath(genpath(data_path));

config_path = fullfile(main_path,'/../config/');
addpath(genpath(config_path));

src_path = fullfile(main_path,'/../src/');
addpath(genpath(src_path));

%output_path = fullfile(main_path,'/../data_output/');
%addpath(genpath(output_path));

georb_path = fullfile(main_path,'/../');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration file for setting all of the required parameters and data
config_filename = 'main_config.in';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call the source code package
georb_function(config_filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cd(georb_path);

% Remove the added paths through restore to the default ones
restoredefaultpath;
