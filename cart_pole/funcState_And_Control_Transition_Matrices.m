function [F_x, F_u] = funcState_And_Control_Transition_Matrices(x, u, du, dt)


global m;
global M;
global l;
global g;

x1 = x(1,1);
x2 = x(2,1);
x3 = x(3,1);
x4 = x(4,1);

u = u(1,1);

% Analytical computation of the required matrices
% F_x = zeros(4,4);
% F_x(1,3) = 1;
% F_x(2,4) = 1;
% F_x(3,2) = m*(-2*g*m + l^2*(m - 4*M)*x4^2*cos(x2) + 2*g*(m + 2*M)*cos(2*x2) - l^2*m*x4^2*cos(3*x2) - 4*u*sin(2*x2))/4/(M + m*sin(x2)^2)^2;
% F_x(3,4) = -2*m*l^2*x4*sin(x2)/(M + m*sin(x2)^2);
% F_x(4,2) = (1/(l*(M + m*sin(x2)^2)^2))*(1/2*(-m - 2*M + m*cos(2*x2))*(-g*(m + M)*cos(x2) + l^2*m*x4^2*cos(2*x2) + u*sin(x2))...
%     - 2*m*cos(x2)*sin(x2)*(g*(m + M)*sin(x2) + cos(x2)*(u - l^2*m*x4^2*sin(x2))));
% F_x(4,4) = -2*m*l*x4*sin(x2)*cos(x2)/(M + m*sin(x2)^2);
% 
% F_u = zeros(4,1);
% F_u(3,1) = 1/(M + m*sin(x2)^2);
% F_u(4,1) = cos(x2)/(l*(M + m*sin(x2)^2));

% Numerical computation of the required matrices using forward difference
% method
F_x = zeros(4,4);
F_x(1,3) = 1;
F_x(2,4) = 1;

[F3_x2, F3_x4, F3_u] = funcFirst_Partial_Derivatives_3var(@funcF3_Dynamics,x2,x4,u);
F_x(3,2) = F3_x2;
F_x(3,4) = F3_x4;

[F4_x2, F4_x4, F4_u] = funcFirst_Partial_Derivatives_3var(@funcF4_Dynamics,x2,x4,u);
F_x(4,2) = F4_x2;
F_x(4,4) = F4_x4;

F_u = zeros(4,1);
F_u(3,1) = F3_u;
F_u(4,1) = F4_u;

end