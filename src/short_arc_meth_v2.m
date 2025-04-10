function [orbit_arc_long, orbit_arc_short, N_short_arcs, orbit_model_struct] = short_arc_meth_v2 (orbit_config_fname, N_short_arcs, iveq, orbit_arc_long, orbit_arc_short, edge_offset, orbit_model_struct)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: short_arc_meth.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Short arcs apporach for the initial orbit parameter estimation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Input confiugration file name *.in in format 
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                 29 August 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 07/04/2025  Thomas Loudis Papanikolaou
%             Source Code minor upgrade 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ORBARC_glb = orbit_model_struct.orbit_arc_length_sec; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N_short_arcs > 0 && iveq == -1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify :: orbit arc length    
    orbit_arc_long = ORBARC_glb - edge_offset;
    % orbital period
    T_orb = 5400;
    short_arc_meth = 0;
    %N_shortarc = 0;
    if orbit_arc_long > T_orb
        short_arc_meth  = 1;
        %N_shortarc      = N_short_arcs;
        orbit_arc_short = T_orb;
    else
        orbit_arc_short = orbit_arc_long; 
        N_short_arcs = 0;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit integration iteration number test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if short_arc_meth == 1
if N_short_arcs > 0    
    if iveq < N_short_arcs
        ORBARC_glb = orbit_arc_short;
    else
        % Restore initial values to parameters
        % Restore :: orbit arc length
        ORBARC_glb = orbit_arc_long;        
    end        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update orbit_model struct matrix
orbit_model_struct.orbit_arc_length_sec = ORBARC_glb; 
