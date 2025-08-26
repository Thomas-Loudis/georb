function [GM,ae,Cnm,Snm,sCnm,sSnm, n_max_eqm, m_max_eqm, n_max_veq, m_max_veq, tide_system, gfm_struct] = prm_grav_model(cfg_fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: prm_gfm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Data reading and preprocessing: Read gravity model file and form
%  spherical harmonic coefficients matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - cfg_fname:          Input confiugration file name *.in  
% 
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Loudis Papanikolaou                                    27 May 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified
% 07/06/2022  Thomas Loudis Papanikolaou
%             Code minor modifications
% 30/10/2022  Dr. Thomas Papanikolaou
%             Read orbit configuration format via structure array or file
% 18/05/2025  TLP
%             Minor changes for optiizing the gravity parameter estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cfg_mode = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read orbit configuration structure array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Gravity Field model
    param_keyword = 'Gravity_field_terms';
    [gravity_field_terms] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'gravity_model_fname';
    [gravity_model_fname] = read_param_cfg(cfg_fname,param_keyword);

    param_keyword = 'gravity_model_degree';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    n_max_gfm = sscanf(param_value,'%d %*');
    
    param_keyword = 'gravity_model_order';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    m_max_gfm = sscanf(param_value,'%d %*');
    
% Variational Equations truncation d/o    
    param_keyword = 'veq_gravity_model_degree';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    VEQ_n_max_gfm = sscanf(param_value,'%d %*');

    param_keyword = 'veq_gravity_model_order';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    VEQ_m_max_gfm = sscanf(param_value,'%d %*');

% Gravity Field parameters estimation
    param_keyword = 'grav_field_paramestim_yn';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    grav_field_paramestim_yn = param_value;

    param_keyword = 'grav_paramestim_degree_min';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    grav_paramestim_degree_min = sscanf(param_value,'%d %*');

    param_keyword = 'grav_paramestim_degree_max';
    [param_value] = read_param_cfg(cfg_fname,param_keyword);
    grav_paramestim_degree_max = sscanf(param_value,'%d %*');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read .in configuration file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg_mode == 1
fid = fopen(cfg_fname);
while (~feof(fid))
    line_ith = fgetl(fid);
    str1 = sscanf(line_ith,'%s %*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Gravity Field model
    test = strcmp(str1,'Gravity_field_terms');
    if test == 1
      gravity_field_terms = sscanf(line_ith,'%*s %s %*') ;
    end

    test = strcmp(str1,'gravity_model_fname');
    if test == 1
      gravity_model_fname = sscanf(line_ith,'%*s %s %*') ;
    end
    
    test = strcmp(str1,'gravity_model_degree');
    if test == 1
      n_max_gfm = sscanf(line_ith,'%*s %d %*');
    end

    test = strcmp(str1,'gravity_model_order');
    if test == 1
      m_max_gfm = sscanf(line_ith,'%*s %d %*');
    end    
    
% Variational Equations truncation d/o    
    test = strcmp(str1,'veq_gravity_model_degree');
    if test == 1
      VEQ_n_max_gfm = sscanf(line_ith,'%*s %d %*');
    end

    test = strcmp(str1,'veq_gravity_model_order');
    if test == 1
      VEQ_m_max_gfm = sscanf(line_ith,'%*s %d %*');
    end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sigma_shc = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static gravity field models
test = strcmp(gravity_field_terms,'static');
if test == 1
    [GM,ae,Cnm,Snm,sCnm,sSnm,nmax,tide_system] = gfc(gravity_model_fname, n_max_gfm, sigma_shc); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-variable gravity field models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epoch's MJD
% Call prm_ic for IC and MJDo
[orbit_arc_length, IC_MJDo, IC_Zo_vec, EOP_data, EOP_interp_no] = prm_ic(cfg_fname);
MJDo = IC_MJDo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Models with time-dependent coefficients e.g. GOCO05s, EIGEN-6S4v2
test = strcmp(gravity_field_terms,'time-variable');
if test == 1        
    % EIGEN models format of Time-Variable gravity models
    gfm_name_start = sscanf(gravity_model_fname,'%5c %*'); 
    test_tv = strcmp(gfm_name_start,'EIGEN');
    if test_tv == 1  
        % EIGEN series format of time-variable gravity models
        [GM,ae,Cnm,Snm,sCnm,sSnm,nmax,tide_system] = gfc_tv2(gravity_model_fname, n_max_gfm, sigma_shc, MJDo);       
    else
        % GOCO series format of time-variable models
        [GM,ae,Cnm,Snm,sCnm,sSnm,nmax,tide_system] = gfc_tv1(gravity_model_fname, n_max_gfm, sigma_shc, MJDo);
%         [sec,day,month,year] = MJD_inv(MJDo);
%         sec00_mean = 0;
%         day15 = 15;
%         [JD_mean,MJD_mean] = MJD_date(sec00_mean,day15,month,year)
%         [GM,ae,Cnm,Snm,sCnm,sSnm,nmax,tide_system] = gfc_tv1(gravity_model_fname, n_max_gfm, sigma_shc, MJD_mean);       
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gravity model truncation degree and order for EQM and VEQ
if n_max_gfm < 0 
    n_max_eqm = nmax;
    m_max_eqm = nmax;
else
    n_max_eqm = n_max_gfm;
    m_max_eqm = m_max_gfm;   
end

n_max_veq = VEQ_n_max_gfm;
m_max_veq = VEQ_m_max_gfm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity Field parameter estimation
grav_paramestim_yn = grav_field_paramestim_yn;
grav_paramestim_degree_range(1,1) = grav_paramestim_degree_min;
grav_paramestim_degree_range(1,2) = grav_paramestim_degree_max;

degree_min = grav_paramestim_degree_min;
degree_max = grav_paramestim_degree_max;
order_min  = 0;
order_max  = degree_max;
S_order_min = 1;

% Initialisation of gravity parameters' coefficients matrices 
[N_param_GRAV, Nparam_C, Nparam_S , C_degree_order, S_degree_order, Cnm_paramestim, Snm_paramestim] = gravity_param_ic(degree_min, degree_max, order_min, order_max, S_order_min);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity signal (case studies only)
gravity_signal_yn = 'n';
test_effect_01 = strcmp(gravity_signal_yn,'y');
if test_effect_01 == 1
    gravity_signal_type = 1;
    if gravity_signal_type == 1 
    % GRAVsimul :: gravity solution simulation
    % gravity_model_fname = 'grav_signal_dgeo.gfc'
    % gravity_model_fname = 'GEORB_Gravity_Solution_DeltaSignal_53351-53357.gfc' 
    % gravity_model_fname = 'GEORB_Gravity_Solution_delta_53351.gfc' 
    gravity_model_fname = 'MAGIC_Level2a_HIS_reference_fields_monthly_mtmshc_HIS_31_20020101_20020131_do_180.gfc'

    grav_signal_degree_max = grav_paramestim_degree_max
    [GM_signal,ae_signal,Cnm_signal,Snm_signal,sCnm_signal,sSnm_signal,nmax_signal,tide_system_signal] = gfc(gravity_model_fname, grav_signal_degree_max, sigma_shc); 

    elseif gravity_signal_type == 2
    % MAGIC simulation
    % gravity_model_fname = 'mtmshc_HIS_20041108_20041208.180' 
    gravity_model_fname = 'shc_L2a_monthly_MAGIC_20041208-20050107_120.gfc'
    % gravity_model_fname = 'shc_L2a_monthly_NGGM_20041208-20050107_120.gfc'
    % gravity_model_fname = 'shc_L2a_monthly_Grace-C-like_20041208-20050107_120.gfc'
    
    grav_signal_degree_max = -1
    [GM_signal,ae_signal,Cnm_signal,Snm_signal,sCnm_signal,sSnm_signal,nmax_signal,tide_system_signal] = gfc(gravity_model_fname, grav_signal_degree_max, sigma_shc); 
    C20_slr = 4.49795080362669e-11;
    C20_magic = Cnm_signal(2+1,0+1);
    Cnm_signal(2+1,0+1) = C20_slr;
    end

    % nmax_signal
    degree_signal = nmax_signal;
    order_signal = nmax_signal;
else
    Cnm_signal = 0;
    Snm_signal = Cnm_signal;
    degree_signal = 0;
    order_signal = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strucutre array
gfm_struct.grav_term  = gravity_field_terms;
gfm_struct.GM         = GM;
gfm_struct.radius     = ae;
gfm_struct.degree     = n_max_gfm;
gfm_struct.degree_eqm = [n_max_eqm m_max_eqm];
gfm_struct.degree_veq = [n_max_veq m_max_veq];
gfm_struct.Cnm        = Cnm;
gfm_struct.Snm        = Snm;
gfm_struct.Cnm_sigma  = sCnm;
gfm_struct.Snm_sigma  = sSnm;
gfm_struct.tide_system = tide_system;

% Gravity Field parameters estimation y/n 
gfm_struct.gravity_solution_yn = grav_paramestim_yn;
% 'param_estim_yn' is set to 'n' (off mode) for the initial orbit determination and 
% during the final step of the POD function it is set to 'y' if gravity_solution_yn is 'y'  
gfm_struct.param_estim_yn = 'n';
% gfm_struct.param_estim_yn = grav_paramestim_yn
gfm_struct.param_estim_degree = grav_paramestim_degree_range;
% C,S coefficients matrices to be estimated as corrections to the apriori gravity model Cnm, Snm
gfm_struct.Cnm_estim  = Cnm_paramestim; 
gfm_struct.Snm_estim  = Snm_paramestim;
% OVerall Number of harmonics coefficients (C,S) to be estimated
gfm_struct.parameters_number = N_param_GRAV;
% Degree and Order values matrix of the C,S coefficients to be estimated
gfm_struct.C_degree_order_estim = C_degree_order;
gfm_struct.S_degree_order_estim = S_degree_order;

% Gravity signal (temporal); Optional, applied only for case studies
gfm_struct.gravity_signal_yn = gravity_signal_yn;
% Signal' delta C,S coefficients matrices
gfm_struct.Cnm_signal  = Cnm_signal; 
gfm_struct.Snm_signal  = Snm_signal;
gfm_struct.degree_signal = [degree_signal order_signal];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

