function [x] = funcSimulate_Dynamics(xo, u_new, Horizon, dt)

global m;
global M;
global l;
global g;

x = xo;

for k = 1:(Horizon-1)

    
    F(1,1) = x(3,k);
    
    F(2,1) = x(4,k);
    
    F(3,1) = (m*g*sin(x(2,k))*cos(x(2,k)) - m*l^2*(x(4,k)^2)*sin(x(2,k)))/(M + m*sin(x(2,k))^2);
    
    F(4,1) = (g/l)*sin(x(2,k)) + F(3,1)*cos(x(2,k))/l;
    %     F(4,1) = ((M + m)*g*sin((M + m*sin(x(2,k)))) - m*l^2*x(4,k)^2*sin(x(2,k))*cos(x(2,k)))/((M + m*sin(x(2,k))^2)*l);
    
    
    G = zeros(4,1);
    
    G(3,1) = 1/(M + m*sin(x(2,k))^2);
    
    G(4,1) = (1/(M + m*sin(x(2,k))^2))*cos(x(2,k))/l;

    x(:,k+1) = x(:,k) + F * dt + G * u_new(:,k) * dt;  
      
end