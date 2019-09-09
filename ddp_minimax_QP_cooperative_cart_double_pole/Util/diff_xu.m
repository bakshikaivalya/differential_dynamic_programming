function [f,fx,fu,fxx,fxu,fuu] = diff_xu(fun,x,u,h,order)

if nargout < 2
   f     = fun(x,u);
else
   n     = size(x,1);

   [f, fz, fzz] = fin_diff(@(z) fun(z(1:n,:,:), z(n+1:end,:,:)),[x;u],h,order);

   fx    = fz (:, 1:n,     :);
   fu    = fz (:, n+1:end, :);
   fxx   = fzz(:, 1:n,     1:n,     :);
   fxu   = fzz(:, 1:n,     n+1:end, :);
   fuu   = fzz(:, n+1:end, n+1:end, :);
   
   if size(f,1) == 1
      fx    = permute(fx,[2 3 1]);
      fu    = permute(fu,[2 3 1]);
      fxx   = permute(fxx,[2 3 4 1]);
      fxu   = permute(fxu,[2 3 4 1]);
      fuu   = permute(fuu,[2 3 4 1]);
   end   
end

