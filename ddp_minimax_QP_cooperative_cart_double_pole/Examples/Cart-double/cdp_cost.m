function l = cdp_cost(E, x, u)
% Instantaneous costs for the cart-pole system

u(isnan(u))  = 0;

% 
l     = E.dt*(E.cu*u.^2 +E.cx(2)*x(2,:).^2 +E.cx(3)*x(3,:).^2 +E.cx(4)*x(4,:).^2 ...
   +E.cx(5)*(x(5,:)-2*pi).^2 +E.cx(6)*(x(6,:)-2*pi).^2);

% sigma = 5;
% l     = E.dt*exp(sigma*(E.cu*u.^2 +E.cx(2)*x(2,:).^2 +E.cx(3)*x(3,:).^2 +E.cx(4)*x(4,:).^2 ...
%    +E.cx(5)*(x(5,:)-2*pi).^2 +E.cx(6)*(x(6,:)-2*pi).^2));


%%%%%---- cost options
% +E.cx(5)*(x(5,:) - 2*pi).^2 +E.cx(6)*x(6,:).^2);
% +E.cx(5)*(cos(x(5,:))).^2 +E.cx(6)*(cos(x(6,:))).^2);
% +E.cx(5)*(x(5,:) - 2*pi).^2 -E.cx(6)*cos(x(6,:)));
% -E.cx(5)*cos(x(5,:)) -E.cx(6)*cos(x(6,:)));
%   +E.cx(5)*x(5,:).^2 +E.cx(6)*x(6,:).^2); 
%    -E.cx(5)*cos(x(5,:)) -E.cx(6)*cos(x(6,:)));   %
% +E.cx(5)*(x(5,:) - 2*pi).^2 -E.cx(6)*(x(6,:) - 2*pi).^2);
%  +E.cx(5)*sin(x(5,:))-E.cx(5)*cos(x(5,:)) +E.cx(6)*sin(x(6,:))-E.cx(6)*cos(x(6,:)));