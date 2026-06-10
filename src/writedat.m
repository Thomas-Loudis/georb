function [formatarray] = writedat(outfilename,wno,prno,dat) 

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
% function : writedat.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Purpose: 
% Write computations to output ascii files according to the selected write 
% parameters (width, digits number,...) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Thomas D. Papanikolaou,  AUTH                                  June  2011 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
% Width Number for 1st column: 
if wno > 10 
    wno1st = 10; 
else 
    wno1st = wno; 
end 
 
% Precision matrix : Number of digits 
[sz1 sz2] = size(prno); 
 
% format argument deformation 
formatarray = ''; 
for i = 1 : sz1 
    if i == 1 
        formati = ['%' num2str(wno1st) 'd'];  
    else 
        formati = ['%' num2str(wno) '.' num2str(prno(i,1)) 'f']; 
    end 
    formatarray = [formatarray formati]; 
end 
formatarray = [formatarray ' \n']; 
%formatarray = [formatarray '\r\n']; 
 
% Write data to file 
fid = fopen(outfilename,'w'); 
fprintf(fid,formatarray,dat); 
fclose(fid); 
