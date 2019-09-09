function  [l0, l_x, l_xx ,l_u ,l_uu ,l_ux] = funcCostfunc_Partial_Derivatives(x, u, k, Q, R, dt)

l0 = 0.5*(u'*R*u) + 0.5*(x'*Q*x);
l_x = Q*x;
l_xx = Q;
l_u = R * u;
l_uu = R;
l_ux = zeros(1,4);

end