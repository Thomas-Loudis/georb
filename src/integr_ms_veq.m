function [orbc,err,veqZarray,veqParray] = integr_ms_veq(zo,arc,MSparam,RKparam,eopdat,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration of Equation of Motion and Variational Equations based on
% Multistep methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments
% - zo:           Initial epoch and state vector in Celestial Reference
%                 System GCRS
%   zo = [to ro' vo']
%   to:           MJD of initial epoch in days
%   ro:           initial position in GCRS (m)
%   vo:           initial velocity in GCRS (m/sec)
% - arc:          Orbit arc lenth in seconds
% - RKparam:      Runge-Kutta & Runge-Kutta-Nystrom methods parameters
%   RKparam(1,1): RK method's ID number
%   RKparam(2,1): Stepsize (h)
%   RKparam(3,1): RKN7(6)-8 parameter for stepsize control (lamda)
%   RKparam(4,1): RKN6(4)-4FD parameter for interpolation points
%                 RKparam(4,1) = sigma, 0<sigma<h
% - MSparam:      Multistep method's parameters
%   MSparam(1,1): Multistep method's ID number
%   MSparam(2,1): Order 
%   MSparam(3,1): Stepsize 
%   MSparam(4,1): Start integrator method for the first epochs
% - eop:          Earth Orientation Parameters (EOP) data that are required
%                 for the orbit arc length
% - dpint:        Number of data points (days) that are required for the 
%                 EOP interpolation to the computation epoch
%
% Output arguments:
% - z_ITRS:       State vector in the Terrestrial Reference System ITRS
%   z_ITRS = [t r_ITRS' v_ITRS']
% - orbc:         State vector in the Celestial Reference System GCRS
%   orbc = [ t r_GCRS' v_GCRS' er' ev']
%   t:            Epoch in seconds of TT (Terrestrial Time)
% - err:          Integrator's local truncation errors  err = [er' ev']
%   er:           local truncation error of position
%   ev:           local truncation error of velocity
% - veqZarray:    VEQ array of state transition matrix (6*Epochs x 6)
% - veqParray:    VEQ array of sensitivity matrix      (6*Epochs x np)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                           May 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 04/06/2012  Replacement of the functions accel.m and pdv_accl.m by the
%             combined function veq_accl.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Output arguments currently set to 0 :
err = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multistep methods Coefficients computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of coefficients gj,bmj for Adams-Bashforth method and
% coefficients g*j and b*mj for Adams-Bashforth-Moulton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Indepedent coefficients gj and g*j are computed with recurence relation.
% Coefficients gj and g*j start from index zero (0).
%  gj = [g1...gj...gm]
% In Matlab the equivalent matrix g starts with index 1.
%  Examble: go = g(1,1)
% Numerical value for go: go = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients gj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical value for go: go = 1
gj(1,1) = 1;
for j = 1 : MSparam(2,1)
    sum_g = 0;
    for k = 0 : j-1
        sum_g = sum_g + (1 / (j+1-k)) * gj(1,k+1);
    end
    gj(1,j+1) = 1 - sum_g;
    clear sum_g
end
clear j k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients bmj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : MSparam(2,1)
    sum_b = 0;
    for l = MSparam(2,1)-j : MSparam(2,1)-1
        sum_b = sum_b + gj(1,l+1) * binom_coeff(l,MSparam(2,1)-j);
    end
    bmj(1,j) = (-1)^(MSparam(2,1)-j) * sum_b;
    clear sum_b
end
clear j l
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients g*j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical value for g*o: g_o = 1
gj_ast(1,1) = 1;
for j = 1 : MSparam(2,1)+1
    sum_g = 0;
    for k = 0 : j-1
        sum_g = sum_g + (1 / (j+1-k)) * gj_ast(1,k+1);
    end
    gj_ast(1,j+1) = - sum_g;
    clear sum_g
end
clear j k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients b*mj (bmj_ast)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : MSparam(2,1)
    sum_b = 0;
    for l = MSparam(2,1)-j : MSparam(2,1)-1
        sum_b = sum_b + gj_ast(1,l+1) * binom_coeff(l,MSparam(2,1)-j);
    end
    bmj_ast(1,j) = (-1)^(MSparam(2,1)-j) * sum_b;
    clear sum_b
end
clear j l
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients delta (dj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 0 : MSparam(2,1)+1
    dj(1,j+1) = (1 - j) * gj_ast(1,j+1);
end
clear j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients delta* (d*j: dj_ast)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  d*j for j=0: d*0 = 1
dj_ast(1,1) = 1;
for j = 1 : MSparam(2,1)+1
    dj_ast(1,j+1) = dj(1,j+1) - dj(1,j);
end
clear j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Integrator for the first steps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of first steps are defined according to multistep method's order.
arcRK = (MSparam(2,1) - 1) * MSparam(3,1);
% RK integrator
RKparam(1,1) = MSparam(4,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Integration for the first steps
[orbcRK,errRK,veqZarrayRK,veqParrayRK] = integr_rk_veq(zo,arcRK,RKparam,eopdat,dpint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function evaluations for the first steps of function "f":  f = [v' a']
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - v:  Velocity vector in GCRS
% - a:  Acceleration vector in GCRS
[RK RK2] = size(orbcRK);
%clear RK2

% Number of parameters of sensitivity matrix
[sz1, sz2] = size(veqParrayRK);
Nparam = sz2 - 1;

% Matrices preallocation
f = zeros(RK,6);
accl = zeros(RK,3);
veqZr_dv2 = zeros(3,6,RK);
veqPr_dv2 = zeros(3,Nparam,RK);

for istart = 1 : RK
    v = orbcRK(istart,5:7)';
    zRK = orbcRK(istart,1:7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perturbations Acceleration and Partial Derivatives 
    [fx,fy,fz,pdv_acc,pdv_acc_param] = veq_accl(zRK,eopdat,dpint);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f(istart,:) = [v' fx fy fz];
    accl(istart,:) = [fx fy fz];
    clear tmjd r v zRK fx fy fz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VEQ Function evaluations for the first steps 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix (3D) of 2nd Derivative of State Transition Matrix  (veqZr_dv2)
    veqZrk = veqZarrayRK( (istart-1)*6+1 : (istart-1)*6+6, 2 : end);
    veqZr_rk = veqZrk(1:3,:);
    veqZr_dv2_iveq = pdv_acc * veqZr_rk;
    veqZr_dv2(:,:,istart) = veqZr_dv2_iveq;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix (3D) of 2nd Derivative of Sensitivity Matrix  (veqPr_dv2)
    veqPrk = veqParrayRK( (istart-1)*6+1 : (istart-1)*6+6, 2 : end);
    veqPr_rk = veqPrk(1:3,:);
    veqPr_dv2_iveq = pdv_acc * veqPr_rk + pdv_acc_param;
    veqPr_dv2(:,:,istart) = veqPr_dv2_iveq;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions for Multistep methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time conversion to seconds (MJD from days to seconds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz1, sz2] = size(orbcRK);
% Initial Epoch (MJD in seconds)
to = orbcRK(sz1,1) * (24*3600);
tmax = to + arc - arcRK;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Position and Velocity vector in GCRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ro: denoted by "rGCRS"
rGCRS = [orbcRK(sz1,2) orbcRK(sz1,3) orbcRK(sz1,4)]';
% vo: denoted by "vGCRS"
vGCRS = [orbcRK(sz1,5) orbcRK(sz1,6) orbcRK(sz1,7)]';
clear sz1 sz2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VEQ Initial arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz1 sz2] = size(veqZarrayRK);
veqZo_ms = veqZarrayRK(sz1-5 : end , 2 : end);
clear sz1 sz2
[sz1 sz2] = size(veqParrayRK);
veqPo_ms = veqParrayRK(sz1-5 : end, 2 : end);
clear sz1 sz2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multistep method is defined by the element "MSparam(1,1)"
% Integration step size
h = MSparam(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices Preallocation
Nepochs = fix(arc/h) + 1;
orbc = zeros(Nepochs,7);
arc_ms = arc - arcRK;
orbc_ms = zeros(fix(arc_ms/h), 7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_orbc_ms = arc_ms/h;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Integration based on Multistep Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gauss-Jackson (Second Sum) methods
if MSparam(1,1) == 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
irun = 1;
for t = to : h : tmax-1
% Integration : EQM
    zoGJ = [rGCRS; vGCRS];
    [Va] = gauss_jackson_int(zoGJ,accl,gj_ast,dj_ast,MSparam);
    [z] = gauss_jackson(Va,gj,dj,MSparam);
    clear zoGJ Va
    
% Integration : VEQ
    % State Transition Matrix
    [V_veqZr_dv2] = gauss_jackson_int_veq(veqZo_ms,veqZr_dv2,gj_ast,dj_ast,MSparam);
    [veqZ] = gauss_jackson_veq(V_veqZr_dv2,gj,dj,MSparam);
    %clear V_veqZr_dv2
    % Sensitivity Matrix
    [V_veqPr_dv2] = gauss_jackson_int_veq(veqPo_ms,veqPr_dv2,gj_ast,dj_ast,MSparam);
    [veqP] = gauss_jackson_veq(V_veqPr_dv2,gj,dj,MSparam);    

% Initial Conditions at at next epoch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % State vector at next epoch
    TT = t + h;
    rGCRS = [z(1,1); z(2,1); z(3,1)];
    vGCRS = [z(4,1); z(5,1); z(6,1)];
    orbc_ms(irun,:) = [TT rGCRS' vGCRS'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function evaluations for next point's first steps
    zmjd = [TT/(24*3600) rGCRS' vGCRS'];
    % Acceleration & Partial derivatives
    [fx,fy,fz,pdv_acc,pdv_acc_param] = veq_accl(zmjd,eopdat,dpint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EQM : Function evaluations for next point's first steps
    [na, n2] = size(accl);
    accel_interm = zeros(na,n2);
    for i1 = 1 : na-1
        for i2 = 1 : n2
            accel_interm(i1,i2) = accl(i1+1,i2);
        end
    end
    accel_interm(na,:) = [fx fy fz];
    accl = accel_interm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VEQ (next epoch)
    t6x1 = TT * [1; 1; 1; 1; 1; 1];
    % State Transition matrix
    veqZarrayMS( (irun-1)*6+1 : (irun-1)*6+6 , : ) = [t6x1 veqZ];
    veqZo_ms = veqZ;
    % Sesnitivity matrix
     veqParrayMS( (irun-1)*6+1 : (irun-1)*6+6 , : ) = [t6x1 veqP];
     veqPo_ms = veqP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VEQ Function evaluations for next point's first steps
    % State Transition matrix
    veqZr = veqZ(1:3,:);
    veqZr_dv2_next_epoch = pdv_acc * veqZr;
    [sz1, sz2, sz3] = size(veqZr_dv2);    
    veqZr_dv2_interm = zeros(sz1, sz2, sz3);    

    %veqZr_dv2_interm(:,:,1:sz3-1) = veqZr_dv2(:,:,2:end);
    for i1 = 1 : sz1
        for i2 = 1 : sz2
            for i3 = 1 : sz3-1
                veqZr_dv2_interm(i1,i2,i3) = veqZr_dv2(i1,i2,i3+1);
            end
        end
    end
    veqZr_dv2_interm(:,:,sz3) = veqZr_dv2_next_epoch;
    veqZr_dv2 = veqZr_dv2_interm;
    
    % 2nd Derivative of Sensitivity Matrix  (veqPr_dv2)
    veqPr = veqP(1:3,:);
    veqPr_dv2_next_epoch = pdv_acc * veqPr + pdv_acc_param;
    %veqPr_dv2(:,:,istart) = veqPr_dv2_iveq;
    [sz1, sz2, sz3] = size(veqPr_dv2);    
    veqPr_dv2_interm = zeros(sz1, sz2, sz3);    

    %veqPr_dv2_interm(:,:,1:sz3-1) = veqPr_dv2(:,:,2:end);
    for i1 = 1 : sz1
        for i2 = 1 : sz2
            for i3 = 1 : sz3-1
                veqPr_dv2_interm(i1,i2,i3) = veqPr_dv2(i1,i2,i3+1);
            end
        end
    end
    veqPr_dv2_interm(:,:,sz3) = veqPr_dv2_next_epoch;
    veqPr_dv2 = veqPr_dv2_interm;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    irun = irun + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gauss-Jackson (Second Sum) methods
% "Predictor-Corrector algorithm"
elseif MSparam(1,1) == 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
irun = 1;
for t = to : h : tmax-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Predictor Step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    zoGJ = [rGCRS; vGCRS];
    [Va] = gauss_jackson_int(zoGJ,accl,gj_ast,dj_ast,MSparam);     
    [z] = gauss_jackson(Va,gj,dj,MSparam);                          
    % Predictor Evaluation Step
    rp = [z(1,1); z(2,1); z(3,1)];
    vp = [z(4,1); z(5,1); z(6,1)];
    TT = t + h;
    % Predictor Step : VEQ
    [V_veqZr_dv2] = gauss_jackson_int_veq(veqZo_ms,veqZr_dv2,gj_ast,dj_ast,MSparam);
    [veqZ] = gauss_jackson_veq(V_veqZr_dv2,gj,dj,MSparam);
    clear V_veqZr_dv2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Corrector Step:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Corrector : Function evaluations for first steps    
    zmjd = [TT/(24*3600) rp' vp'];
    % Acceleration & Partial derivatives
    [fx,fy,fz,pdv_acc] = veq_accl(zmjd,eopdat,dpint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [na n2] = size(accl);
    ap = accl(2:na,:);
    clear na n2
    ap = [ap
          fx fy fz];
    clear rp fx fy fz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Corrector : VEQ
    veqZo_ms = veqZ;
%     veqPo_ms = veqP;
    % Corrector : VEQ Function evaluations for first steps
    %[pdv_acc] = pdv_accl(zmjd,eopdat,dpint);
    veqZr = veqZ(1:3,:);
    veqZr_dv2_next = pdv_acc * veqZr;
    veqZr_dv2 = veqZr_dv2(:,:,2:end);
    veqZr_dv2(:,:,end+1) = veqZr_dv2_next;   
    clear veqZ pdv_acc V_veqZr_dv2_next
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Corrector Step:
    zoGJ = z;
    [Va] = gauss_jackson_int(zoGJ,ap,gj_ast,dj_ast,MSparam);
    [z] = gauss_jackson(Va,gj_ast,dj_ast,MSparam);
    clear zoGJ Va ap
    % Corrector Step : VEQ
    [V_veqZr_dv2] = gauss_jackson_int_veq(veqZo_ms,veqZr_dv2,gj_ast,dj_ast,MSparam);
    [veqZ] = gauss_jackson_veq(V_veqZr_dv2,gj,dj,MSparam);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Next Epoch estimate and function evaluations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % State vector for next epoch TT
    rGCRS = [z(1,1); z(2,1); z(3,1)];
    vGCRS = [z(4,1); z(5,1); z(6,1)];
    orbc_ms(irun,:) = [TT rGCRS' vGCRS'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function evaluations for next point's first steps
    zmjd = [TT/(24*3600) rGCRS' vGCRS'];
    % Acceleration & Partial derivatives
    [fx,fy,fz,pdv_acc] = veq_accl(zmjd,eopdat,dpint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [na n2] = size(accl);
    accl = accl(2:na,:);                                        
    clear na n2
    accl = [accl
             fx fy fz];
    clear fx fy fz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VEQ (next epoch)
    t6x1 = TT * [1; 1; 1; 1; 1; 1];
    veqZarrayMS( (irun-1)*6+1 : (irun-1)*6+6 , : ) = [t6x1 veqZ];
    clear t6x1
    veqZo_ms = veqZ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VEQ Function evaluations for next point's first steps
    %[pdv_acc] = pdv_accl(zmjd,eopdat,dpint);
    veqZr = veqZ(1:3,:);
    veqZr_dv2_next = pdv_acc * veqZr;
    veqZr_dv2(:,:,end) = veqZr_dv2_next;
    clear veqZ pdv_acc V_veqZr_dv2_next
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    irun = irun + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epochs MJD conversion from seconds to days
orbc_ms(:,1) = orbc_ms(:,1) / (24 * 3600);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove columns "er" & "ev" from "orbcRK"
orbcRK_2 = [orbcRK(:,1:7)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit collective matrix
orbc = [orbcRK_2; orbc_ms];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VEQ collective arrays 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Epochs MJD conversion from seconds to days
veqZarrayMS(:,1) = veqZarrayMS(:,1) / (24 * 3600);
veqParrayMS(:,1) = veqParrayMS(:,1) / (24 * 3600);

veqZarray = [ veqZarrayRK
              veqZarrayMS ];
          
veqParray = [ veqParrayRK
              veqParrayMS ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
