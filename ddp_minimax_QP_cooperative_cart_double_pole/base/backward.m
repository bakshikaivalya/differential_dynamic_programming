    %====== STEP 2: backward pass, compute optimal control law and cost-to-go
    backPassDone   = 0;
    while ~backPassDone
        
        t_back   = tic;
        [diverge, Vx, Vxx, l, L, dV] = back_pass(cx,cu,cxx,cxu,cuu,fx,fu,fxx,fxu,fuu,lambda,Op.regType,Op.lims,u);
        back_t   = back_t + toc(t_back);
        
        if diverge
            if verbosity > 2
                fprintf('Cholesky failed at timestep %d.\n',diverge);
            end
            dlambda   = max(dlambda * Op.lambdaFactor, Op.lambdaFactor);
            lambda    = max(lambda * dlambda, Op.lambdaMin);
            if lambda > Op.lambdaMax
                break;
            end
            continue
        end
        backPassDone      = 1;
    end
    % check for termination due to small gradient
    g_norm         = mean(max(abs(l) ./ (abs(u)+1),[],1));
    trace(iter,[1 4 7])  = [iter g_norm nan];
    if g_norm < Op.tolGrad && lambda < 1e-5
        dlambda   = min(dlambda / Op.lambdaFactor, 1/Op.lambdaFactor);
        lambda    = lambda * dlambda * (lambda > Op.lambdaMin);
        trace(iter,[2 8])  = [lambda dlambda];
        if verbosity > 0
            fprintf('\nSUCCESS: gradient norm < tolGrad\n');
        end
%         break;
    end