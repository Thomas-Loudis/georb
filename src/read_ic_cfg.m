function [ic_data_object1,ic_mjd_object1, ic_data_objects, ic_mjd_objects] = read_ic_cfg(ic_filename,id_orbitobject)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  read_ic_cfg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Read IC configuration file for getting the IC epochs per orbit object id 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - ic_filename    : Initial Conditions configuration file name
% - id_orbitobject : Orbiting Object/Satellite ID name
%
% Output arguments:
% - satdata_mjd    : Satellite Data line for the input MJD day 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                             20 August 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of IC epochs of input orbit object
N_object_IC_epochs = 0;
N_ic_epochs = 0;
i_endofheader = 0;
fid = fopen(ic_filename,'r');
while (~feof(fid))
    line_i = fgetl(fid);           
    if i_endofheader == 1 && isempty(line_i) == 0        
    line_keyword = sscanf(line_i,'%s %*');
    N_ic_epochs = N_ic_epochs + 1;
    test_keyword = strcmp(line_keyword,id_orbitobject);
    if test_keyword == 1                 
        N_object_IC_epochs = N_object_IC_epochs + 1;
    end  
    end
    keyword1 = sscanf(line_i,'%s %*');
    endofheader = 'end_of_header';
    test = strcmp(keyword1,endofheader);
    if test == 1
        i_endofheader = 1;
    end    
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation/Initialisation of arrays
Nchar_ID = 700;
blank_ith = blanks(1);

ic_data_objects = repmat(blank_ith, N_ic_epochs, Nchar_ID);
ic_mjd_objects = zeros(N_ic_epochs,2);

ic_data_object1 = repmat(blank_ith, N_object_IC_epochs, Nchar_ID);
ic_mjd_object1 = zeros(N_object_IC_epochs,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(ic_filename,'r');
i_endofheader = 0;
i_epochs = 0;
i_object_epochs = 0;
while (~feof(fid))
    line_i = fgetl(fid);            
    if i_endofheader == 1 && isempty(line_i) == 0
    i_epochs = i_epochs + 1;
    line_keyword = sscanf(line_i,'%s %*');
% GRACE-C     ICRF   GPS   calendar   2019 07 18   0.0     1   3317114.29453071626  252279.97141062448 -6033407.00916356966 -6639.98945932193874 -227.41515344014149 -3673.76634420864957
    Date_format = sscanf(line_i,'%*s %*s %*s %s %*');
    test = strcmp(Date_format,'MJD');
    if test == 1
        MJD_day = sscanf(line_i,'%*s%*s%*s%*s %d %*');
        Sec_00h = sscanf(line_i,'%*s%*s%*s%*s%*s %f %*');
    end
    test = strcmp(Date_format,'calendar');
    if test == 1
        year  = sscanf(line_i,'%*s%*s%*s%*s %d %*');
        month = sscanf(line_i,'%*s%*s%*s%*s%*s %d %*');
        day   = sscanf(line_i,'%*s%*s%*s%*s%*s%*s %d %*');
        sec   = sscanf(line_i,'%*s%*s%*s%*s%*s%*s%*s %f %*');
        [JD_ith, MJD_ith] = MJD_date(sec, day, month, year);
        MJD_day = fix(MJD_ith);
        Sec_00h = sec;
    end

    
    %ic_data_objects(i_epochs)   = line_i;    
    ID_Nchar = length(line_i);
    ic_data_objects(i_epochs, 1:ID_Nchar) = sprintf('%c',line_i); 
    ic_mjd_objects(i_epochs,1) = MJD_day;
    ic_mjd_objects(i_epochs,2) = Sec_00h;
    
    test_keyword = strcmp(line_keyword,id_orbitobject);    
    if test_keyword == 1 
    i_object_epochs = i_object_epochs + 1;
    %ic_data_object1 (i_object_epochs)   = line_i;
    ic_data_object1(i_object_epochs, 1:ID_Nchar) = sprintf('%c',line_i); 
    ic_mjd_object1 (i_object_epochs,1) = MJD_day;
    ic_mjd_object1 (i_object_epochs,2) = Sec_00h;        
    end
    
    end
    keyword1 = sscanf(line_i,'%s %*');
    endofheader = 'end_of_header';
    test = strcmp(keyword1,endofheader);
    if test == 1
        i_endofheader = 1;
    end    
end
fclose(fid);
