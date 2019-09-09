function [diverge_u, diverge_v, Vx, Vxx, lu, Lu, lv, Lv, dV] = back_pass_minimax(cx,cu,cv,cxx,cxu,cxv,cuu,cvv,cuv,cvu,fx,fu,fv,lambda,regType)  

if size(cx,1) == 1
   cx    = shiftdim(cx,  1);
   cu    = shiftdim(cu,  1);
   cv    = shiftdim(cv,  1);
   cxx   = shiftdim(cxx, 1);
   cxu   = shiftdim(cxu, 1);
   cxv   = shiftdim(cxv, 1);
   cuu   = shiftdim(cuu, 1);
   cvv   = shiftdim(cvv, 1);
   cuv   = shiftdim(cuv, 1);
   cvu   = shiftdim(cvu, 1);
end

nx     = size(cx,1);
H      = size(cx,2);
nu     = size(cu,1);
nv     = size(cv,1);

lu     = zeros(nu,H-1);
lv     = zeros(nv,H-1);
Lu     = zeros(nu,nx,H-1);
Lv     = zeros(nv,nx,H-1);
Vx     = zeros(nx,H);
Vxx    = zeros(nx,nx,H);
dV     = [0 0];

Vx(:,H)     = cx(:,H);
Vxx(:,:,H)  = cxx(:,:,H);

diverge_u  = 0;
diverge_v  = 0;

for k = H-1:-1:1
    
   Qu  = cu(:,k)      + fu(:,:,k)'*Vx(:,k+1);
   Qv  = cv(:,k)      + fv(:,:,k)'*Vx(:,k+1);
   Qx  = cx(:,k)      + fx(:,:,k)'*Vx(:,k+1);   
   Qux = cxu(:,:,k)'  + fu(:,:,k)'*Vxx(:,:,k+1)*fx(:,:,k) ; %+ p13mm(Vx(:,k+1),fxu(:,:,:,k));
   Qvx = cxv(:,:,k)'  + fv(:,:,k)'*Vxx(:,:,k+1)*fx(:,:,k) ; %+ p13mm(Vx(:,k+1),fxv(:,:,:,k));
   Quu = cuu(:,:,k)   + fu(:,:,k)'*Vxx(:,:,k+1)*fu(:,:,k) ; %+ p13mm(Vx(:,k+1),fuu(:,:,:,k));
   Qvv = cvv(:,:,k)   + fv(:,:,k)'*Vxx(:,:,k+1)*fv(:,:,k) ; %+ p13mm(Vx(:,k+1),fvv(:,:,:,k));
   Quv = cuv(:,:,k)   + fu(:,:,k)'*Vxx(:,:,k+1)*fv(:,:,k) ; %+ p13mm(Vx(:,k+1),fvu(:,:,:,k));
   Qvu = cvu(:,:,k)   + fv(:,:,k)'*Vxx(:,:,k+1)*fu(:,:,k) ; %+ p13mm(Vx(:,k+1),fuv(:,:,:,k));
   Qxx = cxx(:,:,k)   + fx(:,:,k)'*Vxx(:,:,k+1)*fx(:,:,k) ; %+ p13mm(Vx(:,k+1),fxx(:,:,:,k));
   
   VxxF = (Vxx(:,:,k+1) + lambda*eye(nx)*(regType == 2));
   QuxF = cxu(:,:,k)'   + fu(:,:,k)'*VxxF*fx(:,:,k); % + p13mm(Vx(:,k+1),fxu(:,:,:,k));
   QvxF = cxv(:,:,k)'   + fv(:,:,k)'*VxxF*fx(:,:,k); % + p13mm(Vx(:,k+1),fxv(:,:,:,k));
   QuuF = cuu(:,:,k)    + fu(:,:,k)'*VxxF*fu(:,:,k) + lambda*eye(nu)*(regType == 1); % + p13mm(Vx(:,k+1),fuu(:,:,:,k))
   QvvF = cvv(:,:,k)    + fv(:,:,k)'*VxxF*fv(:,:,k);% + lambda*eye(nv)*(regType == 1); % + p13mm(Vx(:,k+1),fvv(:,:,:,k)) 
   QuvF = cuv(:,:,k)    + fu(:,:,k)'*VxxF*fv(:,:,k) ; %+ p13mm(Vx(:,k+1),fvu(:,:,:,k));
   QvuF = cvu(:,:,k)    + fv(:,:,k)'*VxxF*fu(:,:,k) ; %+ p13mm(Vx(:,k+1),fuv(:,:,:,k));
   
   % cholesky decomposition, check for non-PD
   [Ru,du] = chol(QuuF - (QuvF/QvvF)*QvuF);
   if du ~= 0
      diverge_u  = k;
      return;
   end
   
   [Rv,dv] = chol(QvvF - (QvuF/QuuF)*QuvF);
   if dv == 0
      diverge_v  = k;
      return;
   end
   
   % optimal control law
   lLu = -Ru\(Ru'\[Qu-(QuvF/QvvF)*Qv QuxF-(QuvF/QvvF)*QvxF]);
   luk = lLu(:,1);
   Luk = lLu(:,2:nx+1);
   
   lLv = -Rv\(Rv'\[Qv-(QvuF/QuuF)*Qu QvxF-(QvuF/QuuF)*QuxF]);
   lvk = lLv(:,1);
   Lvk = lLv(:,2:nx+1);
   
   % update cost-to-go approximation
   dV          = dV + [luk'*Qu+lvk'*Qv  .5*(luk'*Quu*luk + lvk'*Qvv*lvk + luk'*Quv*lvk + lvk'*Qvu*luk)];

   Vx(:,k)     = Qx  + Luk'*Quu*luk + Lvk'*Qvv*lvk + Luk'*Qu  + Lvk'*Qv  + Qux'*luk + Qvx'*lvk + Luk'*Quv*lvk + Lvk'*Qvu*luk ;
   Vxx(:,:,k)  = Qxx + Luk'*Quu*Luk + Lvk'*Qvv*Lvk + Luk'*Qux + Lvk'*Qvx + Qux'*Luk + Qvx'*Lvk + Luk'*Quv*Lvk + Lvk'*Qvu*Luk ;
%    Vx(:,k)     = Qx   -  QuxF'*(QuuF\Qu)   -  QvxF'*(QvvF\Qv) ;
%    Vxx(:,:,k)  = Qxx  -  QuxF'*(QuuF\QuxF) -  QvxF'*(QvvF\QvxF) ;
   Vxx(:,:,k)  = .5*(Vxx(:,:,k) + Vxx(:,:,k)');
   
   lu(:,k)      = luk; 
   Lu(:,:,k)    = Luk;
   lv(:,k)      = lvk;
   Lv(:,:,k)    = Lvk;
end

% function c = p13mm(a,b)
% 
% c  = permute(sum(bsxfun(@times,a,b),1), [3 2 1]);