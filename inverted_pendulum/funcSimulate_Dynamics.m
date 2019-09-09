function [x] = funcSimulate_Dynamics(xo,u_new,Horizon,dt)

global m;
global l;
global b;
global I;
global g;

x = xo;

for k = 1:(Horizon-1)

    
    Fx(1,1) = x(2,k);
    Fx(2,1) = -(b/I)*x(2,k) + (g/l)*sin(x(1,k));
    
    G_x = [0; 1];
    
    x(:,k+1) = x(:,k) + Fx * dt + G_x * u_new(:,k) * dt;
    
end