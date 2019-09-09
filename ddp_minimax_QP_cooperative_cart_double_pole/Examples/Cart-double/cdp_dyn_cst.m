function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = cdp_dyn_cst(E,x,u,full_DDP)
% combine quadrotor dynamics and cost
% use helper function finite_difference() to compute derivatives

if nargout == 2
    f = cdp_dyn(E,x,u);
    c = cdp_cost(E,x,u);
else
    % state and control indices
    ix = 1:12;
    iu = 13:16;
    
    % dynamics derivatives
    xu_dyn  = @(xu) cdp_dyn(E,xu(ix,:),xu(iu,:));
    J       = finite_difference(xu_dyn, [x; u]);
    fx      = J(:,ix,:);
    fu      = J(:,iu,:);
    
    % cost first derivatives
    xu_cost = @(xu) cdp_cost(E,xu(ix,:),xu(iu,:));
    J       = squeeze(finite_difference(xu_cost, [x; u]));
    cx      = J(ix,:);
    cu      = J(iu,:);
    
    % cost second derivatives
    xu_Jcst = @(xu) squeeze(finite_difference(xu_cost, xu));
    JJ      = finite_difference(xu_Jcst, [x; u]);
    cxx     = JJ(ix,ix,:);
    cxu     = JJ(ix,iu,:);
    cuu     = JJ(iu,iu,:);
    
    % dynamics second derivatives
    if full_DDP
        xu_Jcst = @(xu) finite_difference(xu_dyn, xu);
        JJ      = finite_difference(xu_Jcst, [x; u]);
        JJ      = reshape(JJ, [4 6 size(J)]);
        JJ      = 0.5*(JJ + permute(JJ,[1 3 2 4]));
        fxx     = JJ(:,ix,ix,:);
        fxu     = JJ(:,ix,iu,:);
        fuu     = JJ(:,iu,iu,:);
    else
        [fxx,fxu,fuu] = deal([]);
    end
    
    [f,c] = deal([]);
end