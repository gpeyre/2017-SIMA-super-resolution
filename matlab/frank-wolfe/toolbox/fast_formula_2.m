function A = fast_formula_2(f)

% fast_formula_2 - fast evaluation of a symbolic function of 4 variables
%
%   v = fast_formula_2(f,X1,X2,Y1,Y2)
%
%   f should be a symbolic function of 4 variables.
%
%   return v=f(X1,X2,Y1,Y2) where xx and yy are usual vectors or matrices.
%   
%   Replaces *,/,^ by vectorized .*,./,.^
%
%   Copyright (c) 2017 Gabriel Peyre

syms xx1 xx2 yy1 yy2;
A = char(f(xx1,xx2,yy1,yy2));
A = strrep(A, '*', '.*');
A = strrep(A, '/', './');
A = strrep(A, '^', '.^');
A = strrep(A, 'xx1', 'X1');
A = strrep(A, 'yy1', 'Y1');
A = strrep(A, 'xx2', 'X2');
A = strrep(A, 'yy2', 'Y2');


end