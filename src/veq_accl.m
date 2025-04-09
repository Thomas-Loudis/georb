function [ax,ay,az,pdv_acc,pdv_acc_param] = veq_accl(z,eop,dpint, orbit_model_struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acceleration of Forces model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Computation of the satellite's acceleration due to gravitational and
%  non-gravitational forces.
%
%  Computation of the partial derivatives of the overall acceleration
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
% - pdv_acc  : Partial derivatives of the overall satellite's
%              acceleration in GCRS  (pdv_acc : 3x3 matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                          June 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 26/10/2022  Thomas Loudis Papanikolaou 
%             Code upgrade and replacement by new function force_eqm_veq.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


EQ_mode = 'VEQ';
[ax,ay,az,pdv_acc,pdv_acc_param] = force_eqm_veq(z,eop,dpint, EQ_mode, orbit_model_struct);

