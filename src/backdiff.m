function [Vf] = backdiff(f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backward Differences
% 
% Purpose:
%   Computation of the backward differences of function f that is to be
%   integrated. Backward differences are used in multistep integration
%   methods.
%
% Input arguments
% - f: Fuction evaluations in one component
%      e.g. f = fx or f = ay
%      f = [f1x....fix....fmx]'
%      m: Method's order
%
% Output arguments
% - Vf:  Backward differences of function evaluation fm
%        Vfm = [Vofm V1fm V2fm .......Vnfm],
%        Vofm = fm
%        n = m-1, m: Method's order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         May 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Nf, n2] = size(f);
% clear n2

Vof = f(Nf,1);

% Nf_0 = Nf
% f_0 = f

% Order m
m = Nf;
% Initialisation
f_j = f;
m_j = Nf;

% Vf = zeros(1,Nf-1);
Vf_loop = zeros(1,Nf-1);
for j = 1 : Nf-1
    % Vnf (n = 1....Order-1)
    % Vnf = zeros((Nf-1),1);
    Vnf = zeros((m_j-1),1);

    % for i = 1 : Nf-1
    for i = 1 : m_j-1
        % Column
        Vnf(i,1) = f_j(i+1,1) - f_j(i,1);
    end

    % row: Store element for Vf 

    % [nVnf n2] = size(Vnf);
    % % clear n2
    % j
    % nVnf_j = nVnf

    % nVnf = Nf-1; 
    % Vf(1,j) = Vnf(nVnf,1);
    Vf_loop(1,j) = Vnf(m_j-1,1);

    % Next order of backwards differences
    % clear f
    % f = Vnf;    
    f_j = Vnf;
    % clear Vnf

    % [Nf, n2] = size(f);
    % % clear n2
    % Nf_j = Nf

    % Nf = Nf-1;
    m_j = m_j-1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vf (collective matrix):
Vf = [Vof Vf_loop];