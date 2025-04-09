function [LLN_hn, LLN_ln, LLN_kn] = loadlovenumbers()


earth_elastic_model = 1;
if earth_elastic_model == 1
% loadlovenumbers_fname = 'load_Love_numbers_from_Gegout.250.txt';
loadlovenumbers_fname = 'Load_Love_numbers_CM_Gegout.txt';

LLN_data = load(loadlovenumbers_fname);
[n_max,sz2] = size(LLN_data);

LLN_hn = zeros(n_max, 1);
LLN_ln = zeros(n_max, 1); 
LLN_kn = zeros(n_max, 1);

LLN_hn(:,1) = LLN_data(:,2);
LLN_ln(:,1) = LLN_data(:,3);
LLN_kn(:,1) = LLN_data(:,4); 

elseif earth_elastic_model == 2 
%
% 
end
