    %====== STEP 4: accept (or not), draw graphics
    if fwdPassDone
        
        % print status
        if verbosity > 1
            fprintf('iter: %-3d  cost: %-9.6g  reduction: %-9.3g  gradient: %-9.3g  log10lam: %3.1f\n', ...
                iter, sum(cost(:)), dcost, g_norm, log10(lambda));
        end
        
        % decrease lambda
        dlambda   = min(dlambda / Op.lambdaFactor, 1/Op.lambdaFactor);
        lambda    = lambda * dlambda * (lambda > Op.lambdaMin);
        
        % accept changes
        u              = unew;
        x              = xnew;
        cost           = costnew;
        flgChange      = 1;
        Op.plotFn(x);
        
        % update trace
        trace(iter,:)  = [iter lambda alpha g_norm dcost z sum(cost(:)) dlambda];
        
        % terminate ?
        if dcost < Op.tolFun
            if verbosity > 0
                fprintf('\nSUCCESS: cost change < tolFun\n');
            end
%             break;
        end
        
    else % no cost improvement
        % increase lambda
        dlambda  = max(dlambda * Op.lambdaFactor, Op.lambdaFactor);
        lambda   = max(lambda * dlambda, Op.lambdaMin);
        
        % print status
        if verbosity > 1
            fprintf('iter: %-3d  REJECTED    expected: %-11.3g    actual: %-11.3g    log10lam: %3.1f\n',...
                iter,expected ,dcost, log10(lambda));
        end
        
        % update trace
        trace(iter,:)  = [iter lambda nan g_norm dcost z sum(cost(:)) dlambda];
        
        % terminate ?
        if lambda > Op.lambdaMax,
            if verbosity > 0
                fprintf('\nEXIT: lambda > lambdaMax\n');
            end
%             break;
        end
    end
%     cost_plot  = cost;  cost_plot(:,end) = E.cstfcn(E,x(:,end),zeros(m,1));
    stop           = graphics_ddp(Op.plot,x,u,cost,L,Vx,Vxx,fx,fxx,fu,fuu,trace(1:iter,:),0);