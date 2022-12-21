function [orbc] = integr_ms(zo,arc,MSparam,RKparam,eopdat,dpint)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration of Equation of Motion based on Multistep methods
%
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
% - z_ITRS: State vector in the Terrestrial Reference System ITRS
%   z_ITRS = [t r_ITRS' v_ITRS']
% - orbc: State vector in the Celestial Reference System GCRS
%   orbc = [ t r_GCRS' v_GCRS' er' ev']
%   t:      Epoch in seconds of TT (Terrestrial Time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou, AUTH                                    May  2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multistep methods Coefficients computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of coefficients gj,bmj for Adams-Bashforth method and
% coefficients g*j and b*mj for Adams-Bashforth-Moulton method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
[orbcRK] = integr_rk(zo,arcRK,RKparam,eopdat,dpint);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function evaluations of the first steps of function "f":  f = [v' a']
% - v:  Velocity vector in GCRS
% - a:  Acceleration vector in GCRS
[RK RK2] = size(orbcRK);
clear RK2
for j = 1 : RK
    v = orbcRK(j,5:7)';
    zRK = orbcRK(j,1:7);
    [fx,fy,fz] = accel(zRK,eopdat,dpint);
    f(j,:) = [v' fx fy fz];
    accl(j,:) = [fx fy fz];
    clear tmjd r v zRK fx fy fz
end
clear j RK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Integration based on Multistep Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time conversion to seconds (MJD from days to seconds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sz1 sz2] = size(orbcRK);
% Initial Epoch (MJD in seconds)
to = orbcRK(sz1,1) * (24*3600);
tmax = to + arc - arcRK;
% Initial Position and Velocity vector in GCRS
% ro: denoted by "rGCRS"
rGCRS = [orbcRK(sz1,2) orbcRK(sz1,3) orbcRK(sz1,4)]';
% vo: denoted by "vGCRS"
vGCRS = [orbcRK(sz1,5) orbcRK(sz1,6) orbcRK(sz1,7)]';
clear sz1 sz2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multistep method is defined by the element "MSparam(1,1)"
% Integration step size
h = MSparam(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices Preallocation
Nepochs = fix(arc/h)+1;
orbc = zeros(Nepochs,7);
arc_ms = arc - arcRK;
orbc_ms = zeros(fix(arc_ms/h),7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adams-Bashforth method
if MSparam(1,1) == 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
for t = to : h : tmax-1
    % Integration step
    yo = [rGCRS; vGCRS];
    [y] = AB(yo,f,MSparam,gj,bmj);
    TT = t + h;
    rGCRS = [y(1,1); y(2,1); y(3,1)];
    vGCRS = [y(4,1); y(5,1); y(6,1)];
    orbc_ms(i,:) = [TT rGCRS' vGCRS'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function evaluations for next point's first steps
    zmjd = [TT/(24*3600) rGCRS' vGCRS'];
    [fx,fy,fz] = accel(zmjd,eopdat,dpint);
    [nf n2] = size(f);
    f = f(2:nf,:);
    f = [f
         vGCRS' fx fy fz];
    clear yo y rGCRS vGCRS zmjd nf n2 fx fy fz
    i = i + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adams-Bashforth-Moulton method
% "Predictor-Corrector algorithm"
elseif MSparam(1,1) == 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
for t = to : h : tmax-1
    % Predictor Step: Adams-Bashforth method
    yo = [rGCRS; vGCRS];
    [y] = AB(yo,f,MSparam,gj,bmj);
    % Predictor Evaluation Step
    rp = [y(1,1); y(2,1); y(3,1)];
    vp = [y(4,1); y(5,1); y(6,1)];
    clear y
    TT = t + h;
    % Function evaluations for Adams-Moulton method's first steps
    zmjd = [TT/(24*3600) rp' vp'];
    [fx,fy,fz] = accel(zmjd,eopdat,dpint);    
    [nf n2] = size(f);
    fp = f(2:nf,:);
    fp = [fp
          vp' fx fy fz];
    % Corrector Step: Adams-Moulton method
    [y] = AM(yo,fp,MSparam,gj_ast,bmj_ast);
    clear rp vp zmjd fx fy fz fp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % State vector
    rGCRS = [y(1,1); y(2,1); y(3,1)];
    vGCRS = [y(4,1); y(5,1); y(6,1)];
    orbc_ms(i,:) = [TT rGCRS' vGCRS'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function evaluations for next point's first steps
    zmjd = [TT/(24*3600) rGCRS' vGCRS'];
    [fx,fy,fz] = accel(zmjd,eopdat,dpint);
    f = f(2:nf,:);
    clear nf n2
    f = [f
         vGCRS' fx fy fz];
    clear TT zmjd fx fy fz
    i = i + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MSparam(1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stoermer method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ATTENTION: Velocity is not computed by Stoermer-Cowell methods
%             Set v = r in equations (e.g. vGRCS = rGCRS)
elseif MSparam(1,1) == 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Two initial position vectors
    [n1 n2] = size(orbcRK)
    roi = [orbcRK(n1-1,2:4); orbcRK(n1,2:4)];
    clear n1 n2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
i = 1;
for t = to : h : tmax-1
    % Integration
    [r] = stoermer_cowell(roi,accl,MSparam,dj);
    % Position vector
    TT = t + h;
    rGCRS = r;
    % set v=r
    vGCRS = rGCRS;
    orbc_ms(i,:) = [TT rGCRS' vGCRS'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Next point's Initial Conditions
    roi = [roi(2,:)
           r'      ];    
    % Function evaluations for next point's first steps
    zmjd = [TT/(24*3600) rGCRS' vGCRS'];
    [fx,fy,fz] = accel(zmjd,eopdat,dpint);
    [nf n2] = size(accl);
    accl = accl(2:nf,:);
    accl = [accl
            fx fy fz];
    clear r TT rGCRS vGCRS zmjd nf n2 fx fy fz
    i = i + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stoermer-Cowell method
% "Predictor-Corrector algorithm"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ATTENTION: Velocity is not computed by Stoermer-Cowell methods
%             Set v = r in equations (e.g. vGRCS = rGCRS)
elseif MSparam(1,1) == 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Two initial position vectors
    [n1 n2] = size(orbcRK);
    roi = [orbcRK(n1-1,2:4); orbcRK(n1,2:4)];
    clear n1 n2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
for t = to : h : tmax-1
    % Predictor Step
    [r] = stoermer_cowell(roi,accl,MSparam,dj);
    % Predictor Evaluation Step
    rp = r;
    vp = r;
    TT = t + h;
    % Function evaluations for Cowell method's first steps
    zmjd = [TT/(24*3600) rp' vp'];
    [fx,fy,fz] = accel(zmjd,eopdat,dpint);
    [na n2] = size(accl);
    ap = accl(2:na,:);
    clear na n2
    ap = [ap
          fx fy fz];
    clear rp vp fx fy fz zmjd
    % Corrector Step
    [r] = stoermer_cowell(roi,ap,MSparam,dj_ast);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Position vector
    rGCRS = r;
    vGCRS = rGCRS;
    orbc_ms(i,:) = [TT rGCRS' vGCRS'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Next point's Initial Conditions
    roi = [roi(2,:)
           r'      ];
    % Function evaluation for next point
    zmjd = [TT/(24*3600) rGCRS' vGCRS'];
    [fx,fy,fz] = accel(zmjd,eopdat,dpint);
    [na n2] = size(accl);
    accl = accl(2:na,:);
    clear na n2
    accl = [accl
             fx fy fz];
    clear fx fy fz zmjd
    i = i + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gauss-Jackson (Second Sum) methods
elseif MSparam(1,1) == 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
for t = to : h : tmax-1
    % Integration
    zoGJ = [rGCRS; vGCRS];
    [Va] = gauss_jackson_int(zoGJ,accl,gj_ast,dj_ast,MSparam);
    [z] = gauss_jackson(Va,gj,dj,MSparam);
    clear zoGJ Va
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % State vector for next epoch
    TT = t + h;
    rGCRS = [z(1,1); z(2,1); z(3,1)];
    vGCRS = [z(4,1); z(5,1); z(6,1)];
    orbc_ms(i,:) = [TT rGCRS' vGCRS'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function evaluations for next point's first steps
    zmjd = [TT/(24*3600) rGCRS' vGCRS'];
    [fx,fy,fz] = accel(zmjd,eopdat,dpint);
    % Functions evaluations for next point's first steps
    [na n2] = size(accl);
    accl = accl(2:na,:);
    clear na n2
    accl = [accl
            fx fy fz];
    clear fx fy fz
    i = i + 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gauss-Jackson (Second Sum) methods
% "Predictor-Corrector algorithm"
elseif MSparam(1,1) == 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
for t = to : h : tmax-1
    % Predictor Step
    zoGJ = [rGCRS; vGCRS];
    [Va] = gauss_jackson_int(zoGJ,accl,gj_ast,dj_ast,MSparam);
    [z] = gauss_jackson(Va,gj,dj,MSparam);
    % Predictor Evaluation Step
    rp = [z(1,1); z(2,1); z(3,1)];
    vp = [z(4,1); z(5,1); z(6,1)];
    TT = t + h;
    % Function evaluations for first steps    
    zmjd = [TT/(24*3600) rp' vp'];
    [fx,fy,fz] = accel(zmjd,eopdat,dpint);
    [na n2] = size(accl);
    ap = accl(2:na,:);
    clear na n2
    ap = [ap
          fx fy fz];
    clear rp fx fy fz
    % Corrector Step:
    zoGJ = z;
    [Va] = gauss_jackson_int(zoGJ,ap,gj_ast,dj_ast,MSparam);
    [z] = gauss_jackson(Va,gj_ast,dj_ast,MSparam);
    clear zoGJ Va ap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % State vector for next epoch TT
    rGCRS = [z(1,1); z(2,1); z(3,1)];
    vGCRS = [z(4,1); z(5,1); z(6,1)];
    orbc_ms(i,:) = [TT rGCRS' vGCRS'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function evaluations for next point's first steps
    zmjd = [TT/(24*3600) rGCRS' vGCRS'];
    [fx,fy,fz] = accel(zmjd,eopdat,dpint);
    [na n2] = size(accl);
    accl = accl(2:na,:);
    clear na n2
    accl = [accl
             fx fy fz];
    clear fx fy fz
    i = i + 1;
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