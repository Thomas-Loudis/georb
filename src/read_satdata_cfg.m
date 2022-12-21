function [satdata_mjd] = read_satdata_cfg(infilename,id1_satellite,id2_MJD)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  read_satdata_cfg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read the sat_data.in configuration file for extracting the data file
%  names for input satellite/object and IC epoch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - infilename     : prm.in input file name (configuration file)
% - id1_satellite  : Satellite ID name
% - id2_MJD        : Input MJD day number
%
% Output arguments:
% - satdata_mjd    : Satellite Data line for the input MJD day 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             19 August 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fid = fopen(infilename,'r');
while (~feof(fid))
    line_i = fgetl(fid);    
    line_keyword = sscanf(line_i,'%s %*');
    test_keyword = strcmp(line_keyword,id1_satellite);
    if test_keyword == 1         
        % GRACE-A   calendar 2009 11 17   GNV1B_2009-11-17_A_01.asc    dynamic     GNV1B_2009-11-17_A_01.asc   dynamic      ACC1B_2009-11-17_A_01.asc      SCA1B_2009-11-17_A_01.asc 
        Date_format = sscanf(line_i,'%*s %s %*');

        test = strcmp(Date_format,'MJD');    
        if test == 1
            MJD_ith = sscanf(line_i,'%*s %*s %d %*');
        end

        test = strcmp(Date_format,'calendar');
        if test == 1
            year  = sscanf(line_i,'%*s %*s %d %*');    
            month = sscanf(line_i,'%*s %*s %*s %d %*');
            day   = sscanf(line_i,'%*s %*s %*s %*s %d %*');
            sec = 0;
            [JD_ith, MJD_ith] = MJD_date(sec, day, month, year);
        end

        if id2_MJD == MJD_ith
            satdata_mjd = line_i;
        end
    end        
end
fclose(fid);
