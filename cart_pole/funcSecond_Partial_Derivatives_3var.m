function [f_xx, f_yy, f_zz, f_xy, f_yz, f_zx] = funcSecond_Partial_Derivatives_3var(f,a,b,c)

% Computation of partial derivatives of a function of three variables using
% forward difference method for derivatives of degree 2

% Inputs:
% Function f is either a named or anonymous function and must be inputted
% with the syntax @f in the arguments to this function
% a,b,c are real or complex numbers

f_xx = (f(a + 2*0.01,b,c) - 2*f(a + 0.01,b,c) + f(a,b,c))/0.01^2;

f_yy = (f(a,b + 2*0.01,c) - 2*f(a,b + 0.01,c) + f(a,b,c))/0.01^2;

f_zz = (f(a,b,c + 2*0.01) - 2*f(a,b,c + 0.01) + f(a,b,c))/0.01^2;

f_xy = (f(a + 0.01,b + 0.01,c) - f(a + 0.01,b - 0.01,c) - f(a - 0.01,b + 0.01,c) + f(a,b,c))/(4*0.01^2);

f_yz = (f(a,b + 0.01,c + 0.01) - f(a,b + 0.01,c - 0.01) - f(a,b - 0.01,c + 0.01) + f(a,b,c))/(4*0.01^2);

f_zx = (f(a + 0.01,b,c + 0.01) - f(a + 0.01,b,c - 0.01) - f(a - 0.01,b,c + 0.01) + f(a,b,c))/(4*0.01^2);

end