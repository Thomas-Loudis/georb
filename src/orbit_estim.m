function [Xmatrix,Amatrix2,Wmatrix2] = orbit_estim(orbc, veqZarray, veqParray, obsorbc, COVobs, COVPform)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: orbit_estim.m (extracted from former mainf_DOD.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Orbit Parameters estimation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                         June 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified:
% 15/04/2021, Thomas Papanikolaou
%             Upgrade and rename of function mainf_DOD.m
% 17/08/2022, Thomas Loudis Papanikolaou
%             Upgrade: matrix decomposition approaches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_itr = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Parameters estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation solution covariance matrix hypothesis:
% 1. Weight matrix set to I (hypothesis P=I) (COVPform=0)
% 2. Weight matrix formed by inverse of Covariance matrix
%    GOCE case:
%    - COVPform = 1 >> Diagonal
%    - COVPform = 2 >> sub-Full Covariance matrix defined by variable COVPform via IN file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit parameter estimation: Initial State vector and Force related parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if Nparam_GLOB > 0   
Texclud = -1;
[Xmatrix,Xmatrix2,Wmatrix, Amatrix, Cx_matrix, Cv_matrix] = estimator_orbit (orbc,veqZarray,veqParray,obsorbc,COVobs,Texclud);
%Cv_matrix = 0;
%Pmatrix = 0;
Xprm_estim(i_itr+1,:) = [i_itr Xmatrix'];
%[Wmatrix2] = mjd2mjdtt(Wmatrix,1);
%[Amatrix2] = mjd2mjdtt(Amatrix,1);
Wmatrix2 = Wmatrix;
Amatrix2 = Amatrix;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1<0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write estimated parameters to output file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outfilename = 'ESTIM_Xmatrix.out';
%fid = fopen(outfilename,'w');
fid = fopen(outfilename,'a');
fprintf(fid,'%s\n', 'Parameters Estimation');
fprintf(fid,'%s\n', '---------------------------------------------------------------------------');
fprintf(fid,'%5d',Xprm_estim(i_itr+1,1));
fprintf(fid,'%31.15e ',Xprm_estim(i_itr+1,2:end));
fprintf(fid,'\n');
fclose(fid);
clear outfilename
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
