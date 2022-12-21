function [doodson_multipliers, delaunay_multipliers] = doodson_number(doodson_number_char)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: doodson_number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose:
%  Compute the integer multipliers of the Delaunay variables and the
%  Doodson arguments for a given Doodson number of tidal wave based on the
%  conventions by Doodson (1921)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% - doodson_argument_number          : Dooson number of the tidal frequency 
%
% Output arguments:
% - doodson_multipliers     : 6 Multipliers of the Doodson variables 
% - delaunay_multipliers    : 5 Multipliers of the Delaunay variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dr. Thomas Loudis Papanikolaou                            1 December 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doodson Argument Number' individual digits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ln_number = length(doodson_number_char);
if ln_number < 7
    doodson_number_char = sprintf('%s%s','0',doodson_number_char);
end

char_digit_1 = doodson_number_char(1);
char_digit_2 = doodson_number_char(2);
char_digit_3 = doodson_number_char(3);

char_digit_4 = doodson_number_char(5);
char_digit_5 = doodson_number_char(6);
char_digit_6 = doodson_number_char(7);

% Individual digits as numbers
digit_1 = sscanf(char_digit_1,'%d%*');
digit_2 = sscanf(char_digit_2,'%d%*');
digit_3 = sscanf(char_digit_3,'%d%*');
digit_4 = sscanf(char_digit_4,'%d%*');
digit_5 = sscanf(char_digit_5,'%d%*');
digit_6 = sscanf(char_digit_6,'%d%*');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Doodson arguments' multipliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = digit_1;
n2 = digit_2 - 5;
n3 = digit_3 - 5;
n4 = digit_4 - 5;
n5 = digit_5 - 5;
n6 = digit_6 - 5;
doodson_multipliers = [n1 n2 n3 n4 n5 n6];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delaunay variables' multipliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N1
delaunay_multipliers(1,1) = n4;
% N2
delaunay_multipliers(1,2) = n6;
% N3
delaunay_multipliers(1,3) = n1 - n2 - n3 - n4 - n6;
% N4
delaunay_multipliers(1,4) = n3 + n6;
% N5
delaunay_multipliers(1,5) = n1 - n2 - n3 - n4 + n5 - n6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

