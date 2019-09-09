function [y,c,yx,yu,yxx,yxu,yuu,cx,cu,cxx,cxu,cuu] = combine(DYN, CST, x, u, t)

final    = isnan(u(1,:));
N        = sum(~final);

if nargout == 2
   c  = CST(x, u, t);
   if N > 0
      y  = DYN(x, u, t);
   else
      y  = 0;
   end
else
   [y,yx,yu,yxx,yxu,yuu]   = DYN(x(:,~final), u(:,~final), t);
   [c,cx,cu,cxx,cxu,cuu]   = CST(x,        u,              t);
end

