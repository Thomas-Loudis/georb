function [formatarray] = writeformat(wno,prno)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function : writedat.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Bulid format construction for write computations to output ascii files
%  according to the selected write parameters (width, digits number,...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas D. Papanikolaou                                       August  2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Width Number for 1st column:
if wno > 10
    wno1st = 10;
else
    wno1st = wno;
end

% Precision matrix : Number of digits
[sz1 sz2] = size(prno);

% format argument deformation
formatarray = '';
for i = 1 : sz1
    formati = ['%' num2str(wno) '.' num2str(prno(i,1)) 'f'];
    formatarray = [formatarray formati];
end
formatarray = [formatarray ' \n'];
%formatarray = [formatarray '\r\n'];
