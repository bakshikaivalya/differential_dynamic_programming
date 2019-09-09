    %====== STEP 1: differentiate dynamics and cost along new trajectory
    if flgChange
        t_diff = tic;
        [~,~,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu]   = DYNCST(x, [u nan(m,1)], 1:N+1);
        diff_t = diff_t + toc(t_diff);
        flgChange   = 0;
    end