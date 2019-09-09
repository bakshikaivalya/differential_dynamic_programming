%%
% -- include some paths
try
    rd = '../../';
    addpath(rd);
    addpath([rd 'examples']);
    addpath([rd 'base']);
    addpath([rd 'util']);
catch
end

full_DDP= false;

% --- cart-double pendulum parameters
E.m1 = 0.5;  % [kg]     mass of cart
E.m2 = 0.5;  % [kg]     mass of 1st pendulum
E.m3 = 0.5;  % [kg]     mass of 2nd pendulum
E.l2 = 0.5;  % [m]      length of 1st pendulum
E.l3 = 0.5;  % [m]      length of 2nd pendulum
E.b  = 0.1;  % [Ns/m]   coefficient of friction between cart and ground
E.g  = 9.8;  % [m/s^2]  acceleration of gravity

E.dt     = 0.02;    % time-step
E.RK     = 0;       % use Runge-Kutta ?
T        = 1.2;     % total time

% --- problem parameters
N        = T/E.dt;                % number of time steps
E.cu     = 1e-2;                  % u-cost coefficient
E.cx     = [0 1 1 1 10 10]';      % state cost
x0       = [0 0 0 0 pi pi]';      % initial state
u0       = randn(1,N);            % initial control sequence

% --- select between iLQG and full DDP
order    = 0;

% --- set up finite differencing
h     = 2^-17;    % finite-difference parameter
n     = 6;
m     = 1;

% --- set up dynamics
F     = @(x,u) cdp_dyn(E,x,u);
Fp    = @(x,u) pour(F,x,u);
DYN   = @(x,u,t) diff_xu(Fp, x, u, h, order);

% --- set up cost
C     = @(x,u) cdp_cost(E,x,u);
Cp    = @(x,u) pour(C,x,u);
CST   = @(x,u,t) diff_xu(Cp, x, u, h, order);

% --- combine into DYNCST
DYNCST   = @(x,u,t) combine(DYN,CST,x,u,t);
% DYNCST  = @(x,u,i) cdp_dyn_cst(E,x,u,full_DDP);

% --- options 
Op      = struct('plot',1,'print',3); 
Op.lims = [-100 100]; 

% --- compile forward pass (uncomment to activate)
% build_f_pass(n,m,N,E)
% Op.F_PASS   = @(varargin) cp_forward_pass(varargin{:},E);

%% ---------------------- user-adjustable parameters ------------------------
defaults = {'lims',           [],...            control limits
            'parallel',       true,...          use parallel line-search?
            'Alpha',          10.^linspace(0,-3,8),... backtracking coefficients
            'tolFun',         1e-7,...          reduction exit criterion
            'tolGrad',        1e-5,...          gradient exit criterion
            'maxIter',        50,...           maximum iterations            
            'lambdaInit',     1,...             initial value for lambda
            'dlambdaInit',    1,...             initial value for dlambda
            'lambdaFactor',   1.6,...           lambda scaling factor
            'lambdaMax',      1e10,...          lambda maximum value
            'lambdaMin',      1e-6,...          below this value lambda = 0
            'regType',        1,...             regularization type 1: q_uu+lambda*eye(); 2: V_xx+lambda*eye()
            'zMin',           0,...             minimal accepted reduction ratio
            'plot',           1,...             0: no;  k>0: every k iters; k<0: every k iters, with derivs window
            'print',          2,...             0: no;  1: final; 2: iter; 3: iter, detailed
            'plotFn',         @(x)0,...         user-defined graphics callback
            'cost',           [],...            initial cost for pre-rolled trajectory            
            };