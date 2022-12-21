function [orb_array,orb_xyz] = grace_orb_itsg(filename)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: grace_kinorb_itsg_read 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read GRACE/GRACE-FO kinematic orbit data provided by ITSG at TU-GRAZ
%  and computed as described by Zehentner and Mayer-Guerr (2013).  
% 
%  Zehentner, N.; Mayer-Guerr, T.; 2013, Kinematic orbits for GRACE and
%  GOCE based on raw GPS observations, Presented at the IAG Scientific
%  Assembly 2013, 1-6 September 2013, Potsdam, Germany.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - filename:   GRACE Kinematic orbit data file name
%
% Output arguments:
% - orb_array:  Kinematic Orbit array including all the provided data per epoch
%   orb_array(i,:) = [MJDgps position_vector(x,y,z) veloctiy_vector(Vx,Vy,Vz) accleration_vector(ax,ay,az)]
%   MJDgps:     MJD in GPS time including fraction of the day
% - orb_xyz:    Kinematic Orbit array including only the position vector to be used as pseudo-observations for dynamic orbit determination 
%   orb_xyz(i,:) = [MJDgps position_vector(x,y,z)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Dr. Thomas Papanikolaou                                
% Written:  12 May 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Open file for reading
fid = fopen(filename);
i = 0;
j = 1;

for i = 1 : 6
    line = fgetl(fid);
end
Nepochs = sscanf(line,'%d %*');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of arrays
orb_array = zeros(Nepochs,10);
orb_xyz =   zeros(Nepochs,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_epoch = 0;
while (~feof(fid))
    line = fgetl(fid);
    i_epoch = i_epoch + 1;
    
    mjd = sscanf(line,'%f %*');
    vector = sscanf(line,'%e');    
    r_vector = vector(2:4,1);
    v_vector = vector(5:7,1);
    a_vector = vector(8:10,1);
    
    pos_x = sscanf(line,'%*f %e %*');
    pos_y = sscanf(line,'%*f %*e %e %*');
    pos_z = sscanf(line,'%*f %*e%*e %e %*');

    vel_x = sscanf(line,'%*f %*e%*e%*e %e %*');
    vel_y = sscanf(line,'%*f %*e%*e%*e %*e %e %*');
    vel_z = sscanf(line,'%*f %*e%*e%*e %*e%*e %e %*');
    
    acc_x = sscanf(line,'%*f %*e%*e%*e %*e%*e%*e %e %*');
    acc_y = sscanf(line,'%*f %*e%*e%*e %*e%*e%*e %*e %e %*');
    acc_z = sscanf(line,'%*f %*e%*e%*e %*e%*e%*e %*e%*e %e %*');

    orb_array(i_epoch,:) = [mjd pos_x pos_y pos_z vel_x vel_y vel_z acc_x acc_y acc_z];
    orb_xyz  (i_epoch,:) = [mjd pos_x pos_y pos_z];
end
fclose(fid);
clear i j fid
