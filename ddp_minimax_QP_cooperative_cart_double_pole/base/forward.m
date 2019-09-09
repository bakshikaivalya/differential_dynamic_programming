    %====== STEP 3: line-search to find new control sequence, trajectory, cost
    fwdPassDone  = 0;
    if backPassDone
        t_fwd = tic;
        if Op.parallel  % parallel line-search
            [xnew,unew,costnew] = forward_pass(x0 ,u, L, x(:,1:N), l, Op.Alpha, DYNCST,Op.lims);
            dcost               = sum(cost(:)) - sum(costnew,2);
            [dcost, w]          = max(dcost);
            alpha               = Op.Alpha(w);
            expected            = -alpha*(dV(1) + alpha*dV(2));
            if expected > 0
                z = dcost/expected;
            else
                z = sign(dcost);
                warning('non-positive expected reduction: should not occur');
            end
            if (z > Op.zMin)
                fwdPassDone = 1;
                costnew     = costnew(:,:,w);
                xnew        = xnew(:,:,w);
                unew        = unew(:,:,w);
            end
        else            % serial backtracking line-search
            for alpha = Op.Alpha
                [xnew,unew,costnew]   = forward_pass(x0 ,u+l*alpha, L, x(:,1:N),[],1,DYNCST,Op.lims);
                dcost    = sum(cost(:)) - sum(costnew(:));
                expected = -alpha*(dV(1) + alpha*dV(2));
                if expected > 0
                    z = dcost/expected;
                else
                    z = sign(dcost);
                    warning('non-positive expected reduction: should not occur');
                end
                if (z > Op.zMin)
                    fwdPassDone = 1;
                    break;
                end
            end
        end
        fwd_t = fwd_t + toc(t_fwd);
    end