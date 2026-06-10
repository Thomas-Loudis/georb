%--------------------------------------------------------------------------
%                                 GEORB                                   %
%                                                                         %
%                                                                         %
%              Gravity and precisE ORBit determination system             %
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
% Dr. Thomas Loudis Papanikolaou                                  July 2022
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Citing GEORB:
%
% Thomas Loudis Papanikolaou (2023). GEORB: Release for precise orbit
% determination of low Earth orbiters and satellite gravity missions,
% Software Impacts, doi: https://doi.org/10.1016/j.simpa.2023.100502    
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% georb_main.m script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  georb_main is the main script file that calls the GEORB source code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 05/07/2022, Thomas Loudis Papanikolaou 
%             Upgrade supporting the release as open source software
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code upgrade of version v.1.9.0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all
clc
format long e
fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Directories Paths: source code, configuration files, input data, output data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
main_path = pwd;

data_path = fullfile(main_path,'..','data');
addpath(genpath(data_path));

config_path = fullfile(main_path,'..','config');
addpath(genpath(config_path));

src_path = fullfile(main_path,'..','src');
addpath(genpath(src_path));

output_path = fullfile(main_path,'..','data_output');
addpath(genpath(output_path));

georb_path = fullfile(main_path,'..');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration file for setting all of the required parameters and data
config_filename = 'main_config.in';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call source code package
georb_intro(config_filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove the added paths through restore to the default ones
restoredefaultpath;
