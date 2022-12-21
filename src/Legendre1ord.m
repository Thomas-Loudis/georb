function [dPnm_norm] = Legendre1ord(phi,nmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First order derivatives of Legendre Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% dPnm(costheta) / dtheta , degree n, order m
%
% Computation of the first order derivatives of the Normalized Associated
% Legendre functions
%
% Evaluation of first order derivatives of Legendre functions is based on
% Recurrence relations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou,  AUTH                                   July 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation
dPnm_norm = zeros(nmax+1,nmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized associated Legendre functions
[Pnm_norm] = Legendre_functions(phi,nmax);

theta = pi/2 - phi;
c = cos(theta);
s = sin(theta) ;

% Initial values
dPoo = 0;
dP10 = - sqrt(3) * s ;
dP11 = sqrt(3) * c ;

%-------------------------------------------------------------------------
% Derivatives for order m=0 and variable degree n
for n = 0: nmax
    if n==0
        dPnm_norm(n+1,0+1) = dPoo ;
    elseif n==1
        dPnm_norm(n+1,0+1) = dP10 ;
    else 
        dP_f1 = sqrt(2*n-1) * ( c * dPnm_norm(n-1+1,0+1) - s * Pnm_norm(n-1+1,0+1) ) ;
        dPnm_norm(n+1,0+1) = (sqrt(2*n+1) / n) * ( dP_f1 - ( (n-1) / sqrt(2*n-3) ) * dPnm_norm(n-2+1,0+1) ) ;
    end
end

% Derivatives for equal degree n and order m    
for n = 1 : nmax
    if n == 1
        dPnm_norm(n+1,n+1) = dP11 ; 
    else   
        dPnm_norm(n+1,n+1) = ( sqrt(2*n+1) / sqrt(2*n) ) * ( s * dPnm_norm(n-1+1,n-1+1) + c * Pnm_norm(n-1+1,n-1+1) ) ;
    end
end

% Derivatives for variable degree n and order m
for n = 2 : nmax
    for m = 1 : n
        if n >= (m+1)
            dP_f2 = c * dPnm_norm(n-1+1,m+1) - s * Pnm_norm(n-1+1,m+1) ;
            dP_f3 = sqrt(  ( (n-1+m)*(n-1-m) ) / (2*n-3)  ) * dPnm_norm(n-2+1,m+1) ;
            dPnm_norm(n+1,m+1) = ( sqrt(2*n+1) / sqrt( (n+m)*(n-m) ) ) * ( sqrt(2*n-1) * dP_f2 -  dP_f3 ) ;
        end
    end
end
