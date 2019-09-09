function [f_x, f_y, f_z] = funcFirst_Partial_Derivatives_3var(f,a,b,c)

% Computation of partial derivatives of a function of three variables using
% forward difference method for derivatives of degree 1

% Inputs:
% Function f is either a named or anonymous function and must be inputted
% with the syntax @f in the arguments to this function
% a,b,c are real or complex numbers

f_x = (f(a + 0.001,b,c) - f(a,b,c))/0.001;

f_y = (f(a,b + 0.001,c) - f(a,b,c))/0.001;

f_z = (f(a,b,c + 0.001) - f(a,b,c))/0.001;

end