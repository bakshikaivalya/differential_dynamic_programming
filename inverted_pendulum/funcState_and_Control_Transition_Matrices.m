function [F_x, F_u] = funcState_And_Control_Transition_Matrices(x, u, du, dt)

global m;
global l;
global b;
global I;
global g;

x1 = x(1,1);
x2 = x(2,1);

F_x = [0 1; (g/l)*cos(x1) -(b/l)];
F_u = [0; 1];

end