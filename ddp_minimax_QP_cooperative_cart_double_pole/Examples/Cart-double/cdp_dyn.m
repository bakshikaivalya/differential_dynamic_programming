function xnew = cdp_dyn(E,x,u)
% cart-double pendulum dynamics model
[n,N] = size(x);
%xnew  = zeros(n,N);
dx    = zeros(n,N);

for i=1:N
    dx(:,i) = E.dt * cdp_f(E,x(:,i),u(:,i));
end

if E.RK
   k2    = dt*cp_f(E,x+0.5*k1,u);
   k3    = dt*cp_f(E,x+0.5*k2,u);
   k4    = dt*cp_f(E,x+k3    ,u);
   dx    =(k1+k4)/6+(k2+k3)/3;  
end

xnew  = x + dx;

function dXdt = cdp_f(E,z,u)  %#eml
%
% *Input arguments:*
%
%	
%   z     state                                                    [6 x 1]
%   u     (optional): force u
%
% *Output arguments:*
%   
%   dXdt    if 3 input arguments:      state derivative wrt time
%           if only 2 input arguments: total mechanical energy
%
%   Note: It is assumed that the state variables are of the following order:
%         x:        [m]     position of cart
%         dx:       [m/s]   velocity of cart
%         dtheta1:  [rad/s] angular velocity of inner pendulum
%         dtheta2:  [rad/s] angular velocity of outer pendulum
%         theta1:   [rad]   angle of inner pendulum
%         theta2:   [rad]   angle of outer pendulum
%
%

% set up the system
m1 = E.m1;
m2 = E.m2;
m3 = E.m3;
l2 = E.l2;
l3 = E.l3;
b  = E.b;
g  = E.g;

% if nargin == 3
  
  A = [2*(m1+m2+m3) -(m2+2*m3)*l2*cos(z(5)) -m3*l3*cos(z(6))
       -(3*m2+6*m3)*cos(z(5)) (2*m2+6*m3)*l2 3*m3*l3*cos(z(5)-z(6))
       -3*cos(z(6)) 3*l2*cos(z(5)-z(6)) 2*l3];
  b = [2*u-2*b*z(2)-(m2+2*m3)*l2*z(3)^2*sin(z(5))-m3*l3*z(4)^2*sin(z(6))
       (3*m2+6*m3)*g*sin(z(5))-3*m3*l3*z(4)^2*sin(z(5)-z(6))
       3*l2*z(3)^2*sin(z(5)-z(6))+3*g*sin(z(6))];
  x = A\b;

  dz = zeros(6,1);
  dz(1) = z(2);
  dz(2) = x(1) + 0.000*randn;
  dz(3) = x(2) + 0.000*randn;
  dz(4) = x(3) + 0.000*randn;
  dz(5) = z(3);
  dz(6) = z(4);
  
  dXdt = dz;
  
% else
%   
%   dz = (m1+m2+m3)*z(2)^2/2+(m2/6+m3/2)*l2^2*z(3)^2+m3*l3^2*z(4)^2/6 ...
%        -(m2/2+m3)*l2*z(2)*z(3)*cos(z(5))-m3*l3*z(2)*z(4)*cos(z(6))/2 ...
%        +m3*l2*l3*z(3)*z(4) *cos(z(5)-z(6))/2+(m2/2+m3)*l2*g*cos(z(5)) ...
%        +m3*l3*g*cos(z(6))/2;