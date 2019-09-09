% =========== Initialization =========== %

% --- initial sizes and controls
n   = size(x0, 1);          % dimension of state vector
m   = size(u0, 1);          % dimension of control vector
N   = size(u0, 2);          % number of state transitions
u   = u0;                   % initial control sequence

% --- proccess options
Op  = setOpts(defaults,Op);

verbosity = Op.print;

switch numel(Op.lims)
    case 0
    case 2*m
        Op.lims = sort(Op.lims,2);
    case 2
        Op.lims = ones(m,1)*sort(Op.lims(:))';
    case m
        Op.lims = Op.lims(:)*[-1 1];
    otherwise
        error('limits are of the wrong size')
end

lambda   = Op.lambdaInit;
dlambda  = Op.dlambdaInit;

%% --- initial trajectory
if size(x0,2) == 1
    diverge = true;
    for alpha = Op.Alpha
        [x,un,cost]  = forward_pass(x0(:,1),alpha*u,[],[],[],1,DYNCST,Op.lims);
        % simplistic divergence test
        if all(abs(x(:)) < 1e8)
            u = un;
            diverge = false;
            break
        end
    end
elseif size(x0,2) == N+1 % already did initial fpass
    x        = x0;
    diverge  = false;
    if isempty(Op.cost)
        error('pre-rolled initial trajectory requires cost')
    else
        cost     = Op.cost;
    end
else
    error('pre-rolled initial trajectory must be of correct length')
end

Op.plotFn(x);

if diverge
    [Vx,Vxx, stop]  = deal(nan);
    L        = zeros(m,n,N);
    cost     = [];
    timing   = [0 0 0 0];
    trace    = [1    lambda nan   nan    nan   nan sum(cost(:)) dlambda];
    if verbosity > 0
        fprintf('\nEXIT: Initial control sequence caused divergence\n');
    end
    return
end

flgChange   = 1;
t_total     = tic;
diff_t      = 0;
back_t      = 0;
fwd_t       = 0;
stop        = 0;
dcost       = 0;
z           = 0;
expected    = 0;
trace       = zeros(min(Op.maxIter,1e6),8);
trace(1,:)  = [1    lambda nan   nan    nan   nan sum(cost(:)) dlambda];
             %[iter lambda alpha g_norm dcost z   sum(cost)    dlambda];
L           = zeros(m,n,N);
if verbosity > 0
    fprintf('\n=========== begin DDP/iLQG ===========\n');
end
graphics_ddp(Op.plot,x,u,cost,L,[],[],[],[],[],[],trace(1,:),1);