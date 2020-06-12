% We will implement an iterative LQR controller to do a swing up
% maneuver. After the swing up is complete, we stabilize with the LQR
% controller.

% The AA203 teaching staff based this code in part on a problem set from 
% Berkeley CS287, and the cartpole visualization code was written by Russ Tedrake

clear all; clc; close all;

rng(210); % for some consistency!  original:(203)

animate = true; % change to visualize animation
noise = false; % change to toggle stochastic disturbance

u_lb = -4; %-4
u_ub = 3; %3

advance_x = @sim_cartpole;
dt = 0.1; % we work with discrete time

%% SCP Setup

% our reference point is x=0, theta=pi, xdot = 0; thetadot=0; this has the pole vertically UP
% our reference control input is 0; this is a fixed point of the system

x_ref = [0; pi; 0; 0];
u_ref = 0;

clip_val = 150;
sigma = 0.005;
if noise
    noise_variance = [0;0;sigma^2;sigma^2];
else
    noise_variance = [0.;0.;0.;0.];
end

Q = eye(4); R = eye(1);

[K_inf, P_inf] = lqr_infinite_horizon_solution(Q, R); 

num_swingup_steps = 75; %75
ep_length = 200; %200

u_old = randn(1,num_swingup_steps); %initialize with random actions

d = @linearize_dynamics;

Qf = diag([10000, 10000, 1000, 1000]);
Q = diag([10,10,2,2]);
R = eye(1);

goal_state=[0; pi; 0; 0];
start_state=[0; 0; 0; 0]; %[0; 0; 0; 0]
%scale_rand = 0.1;
%start_state = -scale_rand + 2*scale_rand*rand(4,1); % randomly start each state between [-1,1]

%% Run SCP
[x_scp,u_scp,residual_history] = scp_solution(advance_x,d, Q, R, Qf, u_lb, u_ub, goal_state, start_state, u_old', num_swingup_steps, dt); % Yours to implement

save('scp_out', 'x_scp', 'u_scp') % save the state and controls obtained by SCP 

%% Plot and Simulate

load("scp_out.mat") % load scp trajectory

x = zeros(4,ep_length+1);
x(:,1) = start_state;

for t=1:ep_length
    if t<num_swingup_steps
        % swing up maneuver is planned for fixed horizon, execute actions
        % from swing up, up until that horizon
        u(:,t) = u_scp(t);
    else 
        % after swing up is over, switch to infinite horizon LQR
        u(:,t) =  K_inf * (x(:,t) - goal_state) + u_ref;
    end
        
    u(:,t) = max(min(u(:,t),clip_val),-clip_val);
    w = noise_variance.*randn(4,1);
    x(:,t+1) = advance_x(x(:,t), u(:,t), dt) + w;
end

figure; 
plot(x');
legend("x", "theta", "x dot", "theta dot");

figure; 
plot(u,'--');
legend("u");

figure;
semilogy(objective_history); grid on
title('Objective Convergence')
xlabel('SCP Iteration')
ylabel('Objective f(x^k)')

if animate
    for t=1:ep_length
        cartpole_draw(t,x(:,t));
    end
end
