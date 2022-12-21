function [ax,ay,az] = accel(z,eop,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleration of Forces model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the satellite's acceleration due to gravitational and
%  non-gravitational forces
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - z:      Epoch's MJD and state vector
%   z = [mjd r' v']
%   mjd:    Epoch's MJD in days of TT (Terrestrial Time)
%   r:      Position in Celestial Reference System (GCRS)
%   v:      Velocity in Celestial Reference System (GCRS)
% - eop:    Earth Orientation Parameters (EOP) data that are required for
%           the orbit arc length
% - dpint:  Number of data points (days) that are required for the EOP
%           interpolation to the computation epoch
%
% Output arguments:
% - ax,ay,az : Satellite's acceleration cartesian components in GCRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                      November 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% ../05/2011  upgrade for including planetary and tides perturbations
% 23/06/2012  New approach to Ocean Tides perturbations with FES2004 model
%             for minimizing computations time. 
%             Use of upgraded functions tides_ocean2.m & tides_fes2004_2.m
% 26/10/2022  Thomas Loudis Papanikolaou 
%             Code upgrade and replacement by new function force_eqm_veq.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EQ_mode = 'EQM';
[ax,ay,az,pdv_acc,pdv_acc_param] = force_eqm_veq(z,eop,dpint, EQ_mode);
