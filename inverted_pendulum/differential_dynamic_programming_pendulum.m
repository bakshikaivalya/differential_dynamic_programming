%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Differential Dynamic Programming for the inverted pendulum example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Nonlinear optimal control using an iterative linear-quadratic approach %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Kaivalya Bakshi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

global m;
global l;
global b;
global I;
global g;

% Inverted Pendulum Dynamics Model Parameters
m = 1; % kg
l = 0.5; % m
b = 0.1;
I = m*l^2; % kg*m^2
g = 9.81; % m/s^2

% Time Horizon 
Horizon = 200;
% Number of Iterations
num_iter = 100;

% Time Discretization step size
dt = 0.01;

% Weight in Final State:
Q_f = zeros(2,2);
Q_f(1,1) = 100;
Q_f(2,2) = 1;

% Weight in the Control and State:
R = 0.01*eye(1,1);
Q = 0.01*eye(2,2);

% Initial State:
xo = zeros(2,1);
xo(1,1) = pi;
xo(2,1) = 0;

% Initial Control:
u_k = zeros(1,Horizon-1);
du_k = zeros(1,Horizon-1);

% Initial trajectory:
x_traj = zeros(2,Horizon);

% Target State: 
p_target(1,1) = 0;
p_target(2,1) = 0;

% Learning Rate:
gamma = 0.5;
 
 
for k = 1:num_iter


    %------------------------------------------------> Linear part of second order approximation of the dynamics
    %------------------------------------------------> Quadratic Approximations of the running cost function
for  j = 1:(Horizon-1)
    
    [l0, l_x, l_xx, l_u, l_uu, l_ux] = funcCostfunc_Partial_Derivatives(x_traj(:,j), u_k(:,j), j, R, dt);
    q0(j) = dt * l0;
    q_k(:,j) = dt * l_x;
    Q_k(:,:,j) = dt * l_xx;
    r_k(:,j) = dt * l_u;
    R_k(:,:,j) = dt * l_uu;
    N_k(:,:,j) = dt * l_ux; 
    
    [dfx, dfu] = funcState_And_Control_Transition_Matrices(x_traj(:,j), u_k(:,j), du_k(:,j), dt);
   
    A(:,:,j) = eye(2,2) + dfx * dt;
    B(:,:,j) = dfu * dt;
    
end


    %------------------------------------------------> Terminal Conditions
    Vxx(:,:,Horizon)= Q_f;
    Vx(:,Horizon) = Q_f * (x_traj(:,Horizon) - p_target); 
    V(Horizon) = 0.5 * (x_traj(:,Horizon) - p_target)' * Q_f * (x_traj(:,Horizon) - p_target);


    %------------------------------------------------> Backpropagation of the Value Function
for j = (Horizon-1):-1:1
   
    X(:,:,j) = [-(g/l) * sin(x_traj(1,j)) * Vx(2,j+1), 0; 0, 0] * dt;
%     X(:,:,j) = zeros(2,2) * dt;
    U(:,:,j) = zeros(1,1) * dt;
    W(:,:,j) = zeros(1,2) * dt;
    H = R_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * B(:,:,j) + U(:,:,j);
    G = N_k(:,:,j) + B(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j) + W(:,:,j);   
    g = r_k(:,j) +  B(:,:,j)' * Vx(:,j+1);
    
   
   inv_H = inv(H);
   L_k(:,:,j)= - inv_H * G;
   l_k (:,j) = - inv_H *g;  
   

   Vxx(:,:,j) = Q_k(:,:,j) + A(:,:,j)' * Vxx(:,:,j+1) * A(:,:,j) + L_k(:,:,j)' * H * L_k(:,:,j) + L_k(:,:,j)' * G + G' * L_k(:,:,j) + X(:,:,j);
   Vx(:,j)= q_k(:,j) +  A(:,:,j)' *  Vx(:,j+1) + L_k(:,:,j)' * g + G' * l_k(:,j) + L_k(:,:,j)'*H * l_k(:,j);
   V(:,j) = q0(j) + V(j+1)   +   0.5 *  l_k (:,j)' * H * l_k (:,j) + l_k (:,j)' * g;
   
end


    %----------------------------------------------> Find the controls
    dx = zeros(2,1);
for i=1:(Horizon-1)
    
    du = l_k(:,i) + L_k(:,:,i) * dx;
    
    % Calculation of increment in dynamics using second order approximation
    % of the dynamics
    Od = [-0.5 * (g/l) * sin(x_traj(1,i)) * dx(1,1)^2; 0] * dt; % computation of second order terms of the approximation of the dynamics
    dx = A(:,:,i) * dx + B(:,:,i) * du + Od; % using second order approximation of the dynamics
    u_new(:,i) = u_k(:,i) + gamma * du;
    
end

    u_k = u_new;


    %---------------------------------------------> Simulation of the Nonlinear System
    [x_traj] = funcSimulate_Dynamics(xo, u_new, Horizon, dt);
    [Cost(:,k)] =  funcCostComputation(x_traj, u_k,p_target, dt, Q_f, Q, R);
    x1(k,:) = x_traj(1,:); 

    fprintf('iLQG Iteration %d,  Current Cost = %e \n',k,Cost(1,k));
 
 
end


time(1)=0;
for i= 2:Horizon
time(i) =time(i-1) + dt;  
end

      
%---------------------------------------------> Plots
figure(1);
subplot(2,2,1)
hold on
plot(time,x_traj(1,:),'linewidth',4);  
plot(time,p_target(1,1)*ones(1,Horizon),'red','linewidth',4)
title('\theta (rad)','fontsize',20); 
xlabel('Time in sec','fontsize',20)
hold off;
grid;
   
   
subplot(2,2,2);
hold on;
plot(time,x_traj(2,:),'linewidth',4); 
plot(time,p_target(2,1)*ones(1,Horizon),'red','linewidth',4)
title('d\theta/dt (rad/s)', 'fontsize',20);
xlabel('Time in sec','fontsize',20)
hold off;
grid;
   
subplot(2,2,3);
hold on
plot(Cost,'linewidth',2); 
xlabel('Iterations','fontsize',20)
title('Cost','fontsize',20);
hold off
grid;