function [trs2crs_mat,dtrs2crs_mat, crs2trs_mat,dcrs2trs_mat, UT1, X_PN, Y_PN, s_PN, xp, yp] = trs2crs(mjd,eop,dpint, orbit_model_struct) 

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
 
 
%-------------------------------------------------------------------------- 
% Function: trs2crs 
%-------------------------------------------------------------------------- 
% Purpose 
%  Implementation of the tranformation between Geocentric Celestial 
%  Reference System (GCRS) and International Terrestrial Reference System 
%  (ITRS) based on the IERS Conventions 2003. 
%  The transformation matrices required for the direct and the inverse 
%  tranformation are computed for the position and velocity vectors. 
%  The derivatives of the transformation matrices are being computed. 
% 
%-------------------------------------------------------------------------- 
% Input arguments: 
% - mjd:    computation epoch's MJD in Terrestrial Time (TT) scale 
% - eop:    Earth Orientation Parameters (EOP) data that are required for 
%           the orbit arc length 
% - dpint:  Number of data points (days) that are required for the EOP 
%           interpolation to the computation epoch 
% 
% Output arguments:
% - trs2crs_mat:  ITRS to GCRS transformation matrix Q(t)*R(t)*W(t) 
% - dtrs2crs_mat: Derivative of ITRS to GCRS transformation matrix
% - crs2trs_mat:  GCRS to ITRS transformation matrix  W(-t) * R(-t) * Q(-t) 
% - dcrs2trs_mat: Derivative of GCRS to ITRS transformation matrix
% 
%-------------------------------------------------------------------------- 
% References: 
%  Petit, G.; Luzum, B. IERS Conventions 2010, IERS Technical Note No. 36; 
%  Verlag des Bundesamts für Kartographie und Geodäsie: Frankfurt am Main, 
%  Germany, 2010. 
%--------------------------------------------------------------------------  
% Thomas D. Papanikolaou, AUTH                                November 2007 
%--------------------------------------------------------------------------  
% Last modified: 
%  10/04/2011  Thomas Papanikolaou 
%              Upgrade and function's rename from IERS_EOP.m to iers_rot.m   
%  19/12/2022, Thomas Loudis Papanikolaou 
%              Code minor upgrade and rename to trs2crs. The transformation 
%              model has considered the IERS Conventions 2010 and its 
%              updates 
% 07/04/2025  Thomas Loudis Papanikolaou 
%             Source Code minor upgrade  
% 30/06/2026  Thomas Loudis Papanikolaou 
%             Precesssion-Nutation model minor changes via new code   
% 10/07/2026  Thomas Loudis Papanikolaou 
%             Source Code minor upgrade  
%-------------------------------------------------------------------------- 
 

% Input epoch MJD in TT time scale
mjd_TT = mjd; 

%-------------------------------------------------------------------------- 
% EOP computation for MJD epoch : 
% EOP corrections cosidering tidal variations; 
% IAU 2006/2000A Precesssion-Nutation model
%-------------------------------------------------------------------------- 
[UT1, X_PN, Y_PN, s_PN, xp, yp] = eop_comp(mjd_TT,eop,dpint, orbit_model_struct); 
%-------------------------------------------------------------------------- 

%-------------------------------------------------------------------------- 
% Transformation matrix computations 
%--------------------------------------------------------------------------
% Transformation matrix from ITRS to GCRS :: trs2crs_mat matrix and dtrs2crs_mat is its derivative w.r.t. time 
% Transformation matrix from GCRS to ITRS :: crs2trs_mat matrix and its derivative dcrs2trs_mat
[trs2crs_mat,dtrs2crs_mat,crs2trs_mat,dcrs2trs_mat] = trs2crs_matrix(mjd_TT, UT1, X_PN, Y_PN, s_PN, xp, yp);
%-------------------------------------------------------------------------- 
