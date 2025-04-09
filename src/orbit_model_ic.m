function [orbit_model_struct] = orbit_model_ic (cfg_fname,orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Orbit modelling :: Initial Conditions based on configuration strucutre array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                     7 May 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions and EOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read configuration file :: IC and EOP
[orbit_arc_length_sec, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(cfg_fname);

% Initial epoch' Modified Julian Date number
orbit_model_struct.IC_MJD = IC_MJDo;
% Initial State Vector in ICRF
orbit_model_struct.IC_CRF = [IC_MJDo IC_Zo_vec'];
% Orbit arc length (sec)
orbit_model_struct.orbit_arc_length_sec = orbit_arc_length_sec;
% Earth Orientation Parameters
orbit_model_struct.EOP_data = EOP_data;
orbit_model_struct.EOP_interp_no = EOP_interp_no;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
