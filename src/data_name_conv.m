function [data_filename] = data_name_conv(mission_name_config, orbiting_object_name, datalevel_name, datarelease_ver, MJD)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function:  data_name_conv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Create file name according to data name format conventions based on data type and date
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - mission_name            : Satellite mission name 
% - orbiting_object_name    : Satellite ID name 
% - MJD                     : MJD Number (integer)
%
% Output arguments:
% - param_value    : Parameter value from the configuration file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                           16 November 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data_filename = 'no_datafile';
mission_name  = 'no_mission';
MJD_day = MJD;

mission_test = 'GRACE_FO_mission';
test_mission_keyword = strcmp(mission_name_config, mission_test);
if test_mission_keyword == 1
    mission_name = 'GRACE-FO';
end
mission_test = 'GRACE_mission';
test_mission_keyword = strcmp(mission_name_config, mission_test);
if test_mission_keyword == 1
    mission_name = 'GRACE';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MJD Day Number to Calendar date
[sec,day_no,month_no,year] = MJD_inv(MJD_day);

% Day number
if day_no < 10
    fname_day = sprintf('%s%d','0',day_no);
else
    fname_day = sprintf('%d',day_no);
end    

% Month number
if month_no < 10
    fname_month = sprintf('%s%d','0',month_no);
else
    fname_month = sprintf('%d',month_no);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE/GRACE-FO Level 1b data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE mission :: Accelerometer data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACC Data GRACE Level 1b data (ACC1B_2009-11-17_A_01.asc)
data_keyword = 'ACC1B';
test_keyword = strcmp(datalevel_name,data_keyword);
if test_keyword == 1
    mission_test = 'GRACE-FO';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        dataformat_suffix = '.txt';
    end
    mission_test = 'GRACE';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        test_datarelease_ver_keyword = strcmp(datarelease_ver, '01');
        if test_datarelease_ver_keyword == 1
            % Release 1
            dataformat_suffix = '.asc';
        else
            dataformat_suffix = '.txt';
        end
    end
    Nchar = length(orbiting_object_name);
    id_letter = orbiting_object_name(Nchar);
    release_ver = datarelease_ver;
    format_ext  = dataformat_suffix;
    % ACC file name considering format name conventions
    data_filename = sprintf('%s%1c%d%1c%s%1c%s%1c%s%1c%s%s', data_keyword,'_',year,'-',fname_month,'-',fname_day,'_',id_letter,'_',release_ver,format_ext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-FO mission :: Accelerometer data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACC Data GRACE Level 1b data (ACC1B_2009-11-17_A_01.asc)
data_keyword = 'ACT1B';
test_keyword = strcmp(datalevel_name,data_keyword);
if test_keyword == 1
    mission_test = 'GRACE-FO';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        dataformat_suffix = '.txt';
    end
    mission_test = 'GRACE';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        test_datarelease_ver_keyword = strcmp(datarelease_ver, '01');
        if test_datarelease_ver_keyword == 1
            % Release 1
            dataformat_suffix = '.asc';
        else
            dataformat_suffix = '.txt';
        end
    end
    Nchar = length(orbiting_object_name);
    id_letter = orbiting_object_name(Nchar);
    release_ver = datarelease_ver;
    format_ext  = dataformat_suffix;
    % ACC file name considering format name conventions
    data_filename = sprintf('%s%1c%d%1c%s%1c%s%1c%s%1c%s%s', data_keyword,'_',year,'-',fname_month,'-',fname_day,'_',id_letter,'_',release_ver,format_ext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE/GRACE-FO missions :: Star Camera data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCA Data GRACE Level 1b data (SCA1B_2009-11-17_A_01.asc)
data_keyword = 'SCA1B';
test_keyword = strcmp(datalevel_name,data_keyword);
if test_keyword == 1
    mission_test = 'GRACE-FO';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        dataformat_suffix = '.txt';
    end
    mission_test = 'GRACE';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        test_datarelease_ver_keyword = strcmp(datarelease_ver, '01');
        if test_datarelease_ver_keyword == 1
            % Release 1
            dataformat_suffix = '.asc';
        else
            dataformat_suffix = '.txt';
        end
    end
    Nchar = length(orbiting_object_name);
    id_letter = orbiting_object_name(Nchar);
    release_ver = datarelease_ver;
    format_ext  = dataformat_suffix;
    % ACC file name considering format name conventions
    data_filename = sprintf('%s%1c%d%1c%s%1c%s%1c%s%1c%s%s', data_keyword,'_',year,'-',fname_month,'-',fname_day,'_',id_letter,'_',release_ver,format_ext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE/GRACE-FO missions :: Orbit GNV1b data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GNV1B Data GRACE Level 1b data (GNV1B_2009-11-17_A_01.asc)
data_keyword = 'GNV1B';
test_keyword = strcmp(datalevel_name,data_keyword);
if test_keyword == 1
    mission_test = 'GRACE-FO';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        dataformat_suffix = '.txt';
    end
    mission_test = 'GRACE';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        test_datarelease_ver_keyword = strcmp(datarelease_ver, '01');
        if test_datarelease_ver_keyword == 1
            % Release 1
            dataformat_suffix = '.asc';
        else
            dataformat_suffix = '.txt';
        end
    end
    Nchar = length(orbiting_object_name);
    id_letter = orbiting_object_name(Nchar);
    release_ver = datarelease_ver;
    format_ext  = dataformat_suffix;
    % ACC file name considering format name conventions
    data_filename = sprintf('%s%1c%d%1c%s%1c%s%1c%s%1c%s%s', data_keyword,'_',year,'-',fname_month,'-',fname_day,'_',id_letter,'_',release_ver,format_ext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atmosphere and Ocean De-Aliasing (AOD) effects   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AOD Data :: GRACE Level 1b data e.g AOD1B_2021-07-17_X_06.asc
data_keyword = 'AOD1B';
test_keyword = strcmp(datalevel_name,data_keyword);
if test_keyword == 1
    dataformat_suffix = '.asc';
    release_ver = datarelease_ver;
    format_ext  = dataformat_suffix;
    id_letter = 'X';
    % AOD file name considering format name conventions
    data_filename = sprintf('%s%1c%d%1c%s%1c%s%1c%s%1c%s%s', data_keyword,'_',year,'-',fname_month,'-',fname_day,'_',id_letter,'_',release_ver,format_ext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE/GRACE-FO missions :: KBR data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GNV1B Data GRACE Level 1b data (KBR1B_2021-07-12_Y_04.txt)
data_keyword = 'KBR1B';
test_keyword = strcmp(datalevel_name,data_keyword);
if test_keyword == 1
    mission_test = 'GRACE-FO';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        dataformat_suffix = '.txt';
        id_letter = 'Y';
    end
    mission_test = 'GRACE';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        test_datarelease_ver_keyword = strcmp(datarelease_ver, '01');
        if test_datarelease_ver_keyword == 1
            % Release 1
            dataformat_suffix = '.asc';
        else
            dataformat_suffix = '.txt';
        end
        id_letter = 'X';
    end
    release_ver = datarelease_ver;
    format_ext  = dataformat_suffix;
    % Data file name considering format name conventions
    data_filename = sprintf('%s%1c%d%1c%s%1c%s%1c%s%1c%s%s', data_keyword,'_',year,'-',fname_month,'-',fname_day,'_',id_letter,'_',release_ver,format_ext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-FO missions :: LRI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LRI1B Data GRACE Level 1b data (LRI1B_2021-07-12_Y_04.txt)
data_keyword = 'LRI1B';
test_keyword = strcmp(datalevel_name,data_keyword);
if test_keyword == 1
    mission_test = 'GRACE-FO';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        dataformat_suffix = '.txt';
        id_letter = 'Y';
    end
    release_ver = datarelease_ver;
    format_ext  = dataformat_suffix;
    % Data file name considering format name conventions
    data_filename = sprintf('%s%1c%d%1c%s%1c%s%1c%s%1c%s%s', data_keyword,'_',year,'-',fname_month,'-',fname_day,'_',id_letter,'_',release_ver,format_ext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITSG TU-Graz: Kinematic Orbit Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACEFO-1_kinematicOrbit_2021-07-12.txt
data_keyword = 'kinematicOrbit';
test_keyword = strcmp(datalevel_name,data_keyword);
if test_keyword == 1
    mission_test = 'GRACE-FO';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        dataformat_suffix = '.txt';
        sat_id_keyword = 'GRACE-C';
        test_keyword_satellite = strcmp(orbiting_object_name,sat_id_keyword);
        if test_keyword_satellite == 1
            satellite_name_data = 'GRACEFO-1';
        end
        sat_id_keyword = 'GRACE-D';
        test_keyword_satellite = strcmp(orbiting_object_name,sat_id_keyword);
        if test_keyword_satellite == 1
            satellite_name_data = 'GRACEFO-2';
        end
    end
    mission_test = 'GRACE';
    test_mission_keyword = strcmp(mission_name, mission_test);
    if test_mission_keyword == 1
        dataformat_suffix = '.txt';
        sat_id_keyword = 'GRACE-A';
        test_keyword_satellite = strcmp(orbiting_object_name,sat_id_keyword);
        if test_keyword_satellite == 1
            satellite_name_data = 'GRACE-1';
        end
        sat_id_keyword = 'GRACE-B';
        test_keyword_satellite = strcmp(orbiting_object_name,sat_id_keyword);
        if test_keyword_satellite == 1
            satellite_name_data = 'GRACE-2';
        end
    end    
    format_ext  = dataformat_suffix;
    % Data file name considering format name conventions
    data_filename = sprintf('%s%1c%s%1c%d%1c%s%1c%s%1c%s', satellite_name_data,'_',data_keyword,'_',year,'-',fname_month,'-',fname_day,format_ext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEORB: Orbit Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRACE-C_59412_orbit_crf.orb
data_keyword = 'GEORB_orbit';
test_keyword = strcmp(datalevel_name,data_keyword);
if test_keyword == 1
    satellite_name_data = orbiting_object_name;
    dataformat_suffix = '.orb';
    format_ext  = dataformat_suffix;
    % Data file name considering format name conventions
    %data_filename = sprintf('%s%1c%d%1c%s%s', satellite_name_data,'_',fix(MJD_day),'_','orbit_crf',format_ext);
    data_filename = sprintf('%s%1c%d%1c%s%s', satellite_name_data,'_',fix(MJD_day),'_','orbit_trf',format_ext);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
