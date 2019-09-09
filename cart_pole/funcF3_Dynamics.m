function [F3] = funcF3_Dynamics(x2,x4,u)

global m;
global M;
global l;
global g;

F3 = (m*g*sin(x2)*cos(x2) - m*l^2*(x4^2)*sin(x2))/(M + m*sin(x2)^2) + u*1/(M + m*sin(x2)^2);
    
end