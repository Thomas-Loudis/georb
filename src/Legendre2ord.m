function [d2Pnm_norm] = Legendre2ord(phi,nmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second-order derivatives of Legendre Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% d^2(Pnm(cos(theta))) / dtheta^2 , degree n, order m
%
% Computation of the second-order derivatives of the Normalized Associated
% Legendre functions
%
% Evaluation of second order derivatives of Legendre functions is based on
% Recurrence relations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou,  AUTH                                   July 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocation
d2Pnm_norm = zeros(nmax+1,nmax+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized associated Legendre functions
[Pnm_norm] = Legendre_functions(phi,nmax) ;

% first-order derivatives of normalized associated Legendre functions
[dPnm_norm] = Legendre1ord(phi,nmax) ;

theta = pi/2 - phi;
c = cos(theta);
s = sin(theta) ;

% Initial values
dPoo = 0;
d2Poo = 0;

dP10 = - sqrt(3) * s ;
d2P10 = - sqrt(3) * c ;

dP11 = sqrt(3) * c ;
d2P11 = -sqrt(3) * s ;

%-------------------------------------------------------------------------
% Second order Derivatives for order m=0 and variable degree n
for n = 0: nmax
    if n==0
        d2Pnm_norm(n+1,0+1) = d2Poo ;
    elseif n==1        
        d2Pnm_norm(n+1,0+1) = d2P10 ;
    else      
        d2P_f1 = sqrt(2*n-1) * ( c * d2Pnm_norm(n-1+1,0+1) - 2 * s * dPnm_norm(n-1+1,0+1) - c * Pnm_norm(n-1+1,0+1) ) ;
        d2Pnm_norm(n+1,0+1) = (sqrt(2*n+1) / n) * ( d2P_f1 - ( (n-1) / sqrt(2*n-3) ) * d2Pnm_norm(n-2+1,0+1) ) ;
    end
end

% Second order Derivatives for equal degree n and order m    
for n = 1 : nmax
    if n == 1
        d2Pnm_norm(n+1,n+1) = d2P11 ;    
    else        
        d2P_sum = s * d2Pnm_norm(n-1+1,n-1+1) + 2 * c * dPnm_norm(n-1+1,n-1+1) - s * Pnm_norm(n-1+1,n-1+1) ;
        d2Pnm_norm(n+1,n+1) = ( sqrt(2*n+1) / sqrt(2*n) ) * d2P_sum  ;
    end
end

% Second order erivatives for variable degree n and order m
for n = 1 : nmax
    for m = 1 : n
        if n >= (m+1)
            d2P_f2 = c * d2Pnm_norm(n-1+1,m+1) - 2 * s * dPnm_norm(n-1+1,m+1) - c * Pnm_norm(n-1+1,m+1) ;
            d2P_f3 = sqrt(  ( (n-1+m)*(n-1-m) ) / (2*n-3)  ) * d2Pnm_norm(n-2+1,m+1) ;
            d2Pnm_norm(n+1,m+1) = ( sqrt(2*n+1) / sqrt( (n+m)*(n-m) ) ) * ( sqrt(2*n-1) * d2P_f2 -  d2P_f3 ) ;
        end
    end
end
