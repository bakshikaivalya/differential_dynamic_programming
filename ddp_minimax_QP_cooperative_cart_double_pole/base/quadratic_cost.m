% quadratic cost function
function l = quadratic_cost(x,u,target,Q,R)

N = size(x,2);   l = zeros(1,N);
for i=1:N
    l(i)  = (x(:,i) - target)'*Q*(x(:,i) - target) + u(:,i)'*R*u(:,i);
end