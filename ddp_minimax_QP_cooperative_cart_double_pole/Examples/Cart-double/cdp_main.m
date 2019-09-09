clear all;
close all;
clc;
%% --- cart-doulble pendulum settings
cdp_settings;
%% --- initialization
Initalization;
%% --- solve 
for iter = 1:Op.maxIter
   if stop
      break;
   end
   %% ====== STEP 1: differentiate dynamics and cost along new trajectory
   linearization;
   %% ====== STEP 2: backward pass, compute optimal control law and cost-to-go
   backward;
   %% ====== STEP 3: line-search to find new control sequence, trajectory, cost
   forward;
   %% ====== STEP 4: accept u / v
   find_opt_u;
end

%% --- draw graphics
if stop
    if verbosity > 0
        fprintf('\nEXIT: Terminated by user\n');
    end
end

if iter == Op.maxIter
    if verbosity > 0
        fprintf('\nEXIT: Maximum iterations reached.\n');
    end
end


if ~isempty(iter)
    total_t = toc(t_total);
    if verbosity > 0
        fprintf(['\n'...
            'iterations:   %-3d\n'...
            'final cost:   %-12.7g\n' ...
            'final grad:   %-12.7g\n' ...
            'final lambda: %-12.7e\n' ...
            'time / iter:  %-5.0f ms\n'...
            'total time:   %-5.2f seconds, of which\n'...
            '  derivs:     %-4.1f%%\n'...
            '  back pass:  %-4.1f%%\n'...
            '  fwd pass:   %-4.1f%%\n'...
            '  other:      %-4.1f%% (graphics etc.)\n'...
            '=========== end DDP/iLQG ===========\n'],...
            iter,sum(cost(:)),g_norm,lambda,1e3*total_t/iter,total_t, [diff_t, back_t, fwd_t, (total_t-diff_t-back_t-fwd_t)]*100/total_t);
    end
    trace    = trace(1:max(trace(:,1)),:);
    timing   = [diff_t back_t fwd_t total_t-diff_t-back_t-fwd_t];
    graphics_ddp(Op.plot,x,u,cost,L,Vx,Vxx,fx,fxx,fu,fuu,trace,2); % draw legend
else
    error('Failure: no iterations completed, something is wrong.')
end

%% --- show cost
% cost_x = cdp_cost(E,x,zeros(nx,H+1));
% figure; plot(cost_x);
% title('Final cost');
%% --- play animation
% figure;
% u(:,H+1)  = zeros(nu,1);
% for i = 1:H
%     cdp_draw(x(1,i),x(5,i),x(6,i),u(:,i));
%     pause(E.dt);
% end