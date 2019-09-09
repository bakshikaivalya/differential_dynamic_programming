function [f,fx,fxx] = fin_diff(fun, x, h, order)

if nargout == 1
   f  = fun(x);
   return
end

if nargin < 3
   h  = 2^-17;
end
if nargin < 4
   order = 4;
end

[n N]       = size(x);
p           = (n^2-n)/2;
I           = eye(n);

if nargout == 2
   Q        = [zeros(n,1) I -I];
   D        = [I -I]';
   dX       = h*Q;
   X        = pp(p23(x), dX);
   F        = fun(X);
   f        = F(:,1,:);
   dF       = pp(F(:,2:end,:), -f);
   fx       = mm(dF, D)*(2*h)^-1;
   return
end

T        = triu(ones(n),1);
[n1,n2]  = find(T);
Q1       = I(:,n1);
Q2       = I(:,n2);

switch order
   case 0
      Q        = [zeros(n,1) I -I];
      D        = [[I; -I] [I; I]];
      hh       = [2*h*ones(1,n) h^2*ones(1,n)];
   case 1
      Q        = [zeros(n,1) I -I Q1+Q2];
      D        = [[I I; -I I; zeros(p,2*n)] [-Q1-Q2; zeros(n,p); eye(p)]];
      hh       = [2*h*ones(1,n) h^2*ones(1,n) h^2*ones(1,p)];
   case 2      
      Q        = [zeros(n,1) I -I Q1+Q2 -Q1-Q2];
      D        = [[I I; -I I; zeros(2*p,2*n)] [-Q1-Q2; -Q1-Q2; eye(p); eye(p)]];
      hh       = [2*h*ones(1,n) h^2*ones(1,n) 2*h^2*ones(1,p)];
   otherwise
      Q        = [zeros(n,1) I -I Q1+Q2 -Q1-Q2 Q1-Q2 Q2-Q1];
      B        = (n+1)*I - 2;
      D        = [[I;-I; zeros(4*p,n)]...      
                  [B; B; repmat((Q1+Q2)',[4 1])]...
                  [zeros(2*n, p); eye(p); eye(p); -eye(p); -eye(p)]];
      hh       = [2*h*ones(1,n) (3*n-3)*h^2*ones(1,n) 4*h^2*ones(1,p)];                         
end

dX       = h*Q;
X        = pp(permute(x,[1 3 2]), dX);
F        = fun(X);
m        = size(F,1);
f        = F(:,1,:);
dF       = pp(F(:,2:end,:), -f);
FX       = tt(mms(dF, D), hh.^-1);
fx       = FX(:,1:n,:);
yxx      = FX(:,n+1:end,:);

if order ~= 0
   T(logical(T))        = n+(1:p);
   T                    = T+T';
   T(logical(eye(n)))   = 1:n;
   fxx                  = yxx(:,T(:),:);
else
   T                    = logical(eye(n));
   fxx                  = zeros(m,n*n,N);
   fxx(:,T(:),:)        = yxx;
end

fxx   = reshape(fxx,[m n n N]);

function ab = pp(a,b) 

ab = bsxfun(@plus,a,b);

function ab = tt(a,b)

ab = bsxfun(@times,a,b);
