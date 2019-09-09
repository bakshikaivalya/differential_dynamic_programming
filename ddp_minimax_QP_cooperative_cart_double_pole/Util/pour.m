function y = pour(fun,x,u)

[n N K] = size(x);
m       = size(u,1);
if K > 1
   x = reshape(x, [n N*K]);
   u = reshape(u, [m N*K]);
end

y  = fun(x,u);

y = reshape(y, [size(y,1) N K]);
