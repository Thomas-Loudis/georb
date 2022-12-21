function [config_file] = write_config_struct2file(config_struct,config_file)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  write_config_struct2file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Write configuration structure to output file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - config_struct     : Configuration structure name 
% - config_file     : Configuration file name 
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Papanikolaou                                   5 November 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
% Open file for writing
fid = fopen(config_file,'w');
    
[k m] = size(config_struct);
for i = 1 : m
    % Read Structure
    param_name  = config_struct(1,i).names;
    param_value = config_struct(1,i).values;
    % Write to file
    fprintf(fid,'%-27s %s ',param_name, param_value);
    fprintf(fid,'%s\n','');    
end   
fclose(fid);
