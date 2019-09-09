%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Iterative Linear Quadratic Regulator for Cart Pole Dynamics       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Nonlinear optimal control using an iterative linear quadratic approach %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Kaivalya Bakshi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;


global m;
global M;
global l;
global g;


% Cart Pole Dynamics Model Paramters
m = 1; % kg
M = 10; % kg
l = 0.5; % m
g = 9.81; % m/s^2

% Horizon 
Horizon = 200; % 1.5sec
% Number of Iterations
num_iter = 100;

% Discretization
dt = 0.01;

% Weight in Final State:
Q_f = zeros(4,4);
Q_f(1,1) = 100;
Q_f(2,2) = 100;
Q_f(3,3) = 50;
Q_f(4,4) = 50;

% Weight in the Control and State:
R = 0.01 * eye(1,1);
Q = 0.01  * eye(4,4);

% Initial Configuration:
xo = zeros(4,1);
xo(1,1) = 0;
xo(2,1) = pi;

% Initial Control:
u_k = zeros(1,Horizon-1);
du_k = zeros(1,Horizon-1);


% Initial trajectory:
x_traj = zeros(4,Horizon);
 

% Target: 
p_target(1,1) = 0;
p_target(2,1) = 0;
p_target(3,1) = 0;
p_target(4,1) = 0;


% Learning Rate:
gamma = 0.5;
 
 
for k = 1:num_iter

%------------------------------------------------> Linearization of the dynamics
%------------------------------------------------> Quadratic Approximations of the cost function 
for  j = 1:(Horizon-1)
    
    [l0,l_x,l_xx,l_u,l_uu,l_ux] = funcCostfunc_Partial_Derivatives(x_traj(:,j),u_k(:,j),j,Q,R,dt);
    q0(j) = dt * l0;
    q_k(:,j) = dt * l_x;
    Q_k(:,:,j) = dt * l_xx;
    r_k(:,j) = dt * l_u;
    R_k(:,:,j) = dt * l_uu;
    N_k(:,:,j) = dt * l_ux;
    
    [dfx,dfu] = funcState_And_Control_Transition_Matrices(x_traj(:,j),u_k(:,j),du_k(:,j),dt);
   
    A(:,:,j) = eye(4,4) + dfx * dt;
    B(:,:,j) = dfu * dt;
    
end

%------------------------------------------------> Terminal Conditions
Vxx(:,:,Horizon)= Q_f;
Vx(:,Horizon) = Q_f * (x_traj(:,Horizon) - p_target); 
V(Horizon) = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target);


%------------------------------------------------> Backpropagation of the Value Function
for j = (Horizon-1):-1:1
   
    % Computation of X, U, W
    F1_xx = zeros(4,4);
    F1_uu = zeros(1,1);
    F1_ux = zeros(1,4);
    
    F2_xx = zeros(4,4);
    F2_uu = zeros(1,1);
    F2_ux = zeros(1,4);

    [F3_x2x2, F3_x4x4, F3_uu, F3_x2x4, F3_x4u, F3_ux2] = funcSecond_Partial_Derivatives_3var(@funcF3_Dynamics,x_traj(2,j),x_traj(4,j),u_k(:,j));
    F3_xx = zeros(4,4);
    F3_xx(2,2) = F3_x2x2;
    F3_xx(2,4) = F3_x2x4;
    F3_xx(4,2) = F3_x2x4;
    F3_xx(4,4) = F3_x4x4;
    
    F3_ux = zeros(1,4);
    F3_ux(1,2) = F3_ux2;
    
    F3_uu = zeros(1,1); 
    
   
    [F4_x2x2, F4_x4x4, F4_uu, F4_x2x4, F4_x4u, F4_ux2] = funcSecond_Partial_Derivatives_3var(@funcF4_Dynamics,x_traj(2,j),x_traj(4,j),u_k(:,j));   
    F4_xx = zeros(4,4);
    F4_xx(2,2) = F4_x2x2;
    F4_xx(2,4) = F4_x2x4;
    F4_xx(4,2) = F4_x2x4;
    F4_xx(4,4) = F4_x4x4;
    
    F4_ux = zeros(1,4);
    F4_ux(1,2) = F4_ux2;
        
    F4_uu = zeros(1,1);
    
    
    X(:,:,j) = (F1_xx*Vx(1,j+1) + F2_xx*Vx(2,j+1) + F3_xx*Vx(3,j+1) + F4_xx*Vx(4,j+1)) * dt;
    U(:,:,j) = (F1_uu*Vx(1,j+1) + F2_uu*Vx(2,j+1) + F3_uu*Vx(3,j+1) + F4_uu*Vx(4,j+1)) * dt;
    W(:,:,j) = (F1_ux*Vx(1,j+1) + F2_ux*Vx(2,j+1) + F3_ux*Vx(3,j+1) + F4_ux*Vx(4,j+1)) * dt;

    H = R_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * B(:,:,j);% + U(:,:,j);
    G = N_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j);% + W(:,:,j);  
    g = r_k(:,j) +  B(:,:,j)' * Vx(:,j+1);
   
   
    inv_H = inv(H);
    L_k(:,:,j)= - inv_H * G;
    l_k (:,j) = - inv_H * g;
   
 
   Vxx(:,:,j) = Q_k(:,:,j)+ A(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j) + L_k(:,:,j)' * H * L_k(:,:,j) + L_k(:,:,j)' * G + G' * L_k(:,:,j);% + X(:,:,j);
   Vx(:,j)= q_k(:,j) +  A(:,:,j)' *  Vx(:,j+1) + L_k(:,:,j)' * g + G' * l_k(:,j) + L_k(:,:,j)' * H * l_k(:,j);
   V(:,j) = q0(j) + V(j+1)   +   0.5 *  l_k (:,j)' * H * l_k (:,j) + l_k (:,j)' * g;
end


%----------------------------------------------> Find the controls
dx = zeros(4,1);
for i=1:(Horizon-1)
    
    du = l_k(:,i) + L_k(:,:,i) * dx;
   
   % Calculation of increment in dynamics using second order approximation
   % of the dynamics
   
    % Computation of Od
    F1_xx = zeros(4,4);
    F1_uu = zeros(1,1);
    F1_ux = zeros(1,4);
    
    F2_xx = zeros(4,4);
    F2_uu = zeros(1,1);
    F2_ux = zeros(1,4);

    [F3_x2x2, F3_x4x4, F3_uu, F3_x2x4, F3_x4u, F3_ux2] = funcSecond_Partial_Derivatives_3var(@funcF3_Dynamics,x_traj(2,i),x_traj(4,i),u_k(:,i));
    
    F3_xx = zeros(4,4);
    F3_xx(2,2) = F3_x2x2;
    F3_xx(2,4) = F3_x2x4;
    F3_xx(4,2) = F3_x2x4;
    F3_xx(4,4) = F3_x4x4;
    
    F3_ux = zeros(1,4);
    F3_ux(1,2) = F3_ux2;
    
    F3_uu = zeros(1,1); 
   
    [F4_x2x2, F4_x4x4, F4_uu, F4_x2x4, F4_x4u, F4_ux2] = funcSecond_Partial_Derivatives_3var(@funcF4_Dynamics,x_traj(2,i),x_traj(4,i),u_k(:,i));
    
    F4_xx = zeros(4,4);    
    F4_xx(2,2) = F4_x2x2;
    F4_xx(2,4) = F4_x2x4;
    F4_xx(4,2) = F4_x2x4;
    F4_xx(4,4) = F4_x4x4;
    
    F4_ux = zeros(1,4);
    F4_ux(1,2) = F4_ux2;
    
    F4_uu = zeros(1,1);
    
    
    Od = (0.5*dx'*(F1_xx + F2_xx + F3_xx + F4_xx)*dx + 0.5*du'*(F1_uu + F2_uu+ F3_uu + F4_uu)*du + du'*(F1_ux + F2_ux+ F3_ux + F4_ux)*dx) * dt; % computation of second order terms of the approximation of the dynamics   
         
    dx = A(:,:,i) * dx + B(:,:,i) * du;% + Od;
    
    u_new(:,i) = u_k(:,i) + gamma * du;
    
end

u_k = u_new;


%---------------------------------------------> Simulation of the Nonlinear System
[x_traj] = funcSimulate_Dynamics(xo,u_new,Horizon,dt);
[Cost(:,k)] =  funcCostComputation(x_traj,u_k,p_target,dt,Q_f,Q,R);
x1(k,:) = x_traj(1,:);
 

fprintf('iLQG Iteration %d,  Current Cost = %e \n',k,Cost(1,k));
 
 
end



   time(1)=0;
   for i= 2:Horizon
    time(i) =time(i-1) + dt;  
   end

      
%---------------------------------------------> Plot Section
   figure(1);
   subplot(3,2,1)
   hold on
   plot(time,x_traj(1,:),'linewidth',4);  
   plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4)
   title('x (m)','fontsize',20); 
   xlabel('Time in sec','fontsize',20)
   hold off;
   grid;
   
   
   subplot(3,2,2);hold on;
   plot(time,x_traj(2,:),'linewidth',4); 
   plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4)
   title('\theta (rad)','fontsize',20);
   xlabel('Time in sec','fontsize',20)
   hold off;
   grid;   

    
   subplot(3,2,3);hold on
   plot(time,x_traj(3,:),'linewidth',4); 
   plot(time,p_target(3,1)*ones(1,Horizon),'red','linewidth',4)
   title('dx/dt (m/s)','fontsize',20)
   xlabel('Time in sec','fontsize',20)
   hold off;
   grid;
   
   subplot(3,2,4);hold on
   plot(time,x_traj(4,:),'linewidth',4); 
   plot(time,p_target(4,1)*ones(1,Horizon),'red','linewidth',4)
   title('d \theta/dt (rad/s)','fontsize',20)
   xlabel('Time in sec','fontsize',20)
   hold off;
   grid;
   
   subplot(3,2,5);hold on
   plot(Cost,'linewidth',2); 
   xlabel('Iterations','fontsize',20)
   title('Cost','fontsize',20);