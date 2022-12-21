function [fid] = write_cfg_change(infilename,infilename2,param_keyword,param_write)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  write_cfg_change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Change parameters of input files *.IN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - infilename  : prm.in input file name
% - infilename2 : modified *.in input filename that changes are written 
%
% Output arguments:
% -     : 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                        June  2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 16/10/2022,   Dr. Thomas Loudis Papanikolaou
%               Function revised and renamed 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% fid = fopen(infilename,'r+');
fid = fopen(infilename,'r');
fid2 = fopen(infilename2,'w');
while (~feof(fid))
    lineith = fgetl(fid);
    % Test of keyword
    line_keyword = sscanf(lineith,'%s %*');
    test_keyword = strcmp(line_keyword,param_keyword);
    if test_keyword == 1 
        param_value = sscanf(lineith,'%*s %s %*');
        param_line = sscanf(lineith,'%*s %1000c'); 
        % Update configuration line
        lineith = sprintf('%s ', param_keyword, param_write);
    end  
    
%   fprintf(fid, '%s\n',lineith);
    fprintf(fid2,'%s\n',lineith);
end

fclose(fid);
fclose(fid2);
