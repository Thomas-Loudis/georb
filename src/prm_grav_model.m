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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
end


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
