% initiation of model parameters


clear params;

% size of excavator

% [m]
params.l1 = 10; % length from A to C
params.l2 = 5; %length from C to D
params.l3 = params.l1 + params.l2;
params.l4 = 10; %length of movable excavator arm
params.l5 = 5; %length from shovel to attachment point of blue pulley


params.ang_base = 1/4*pi;
params.g = 9.81;

% [m] cable reel properties
% radius
params.r_B1 = 0.20;
params.r_B2 = 0.20;
params.r_P1 = 0.20;
params.r_P2 = 0.20;

% thickness
params.h_B1 = 0.15;
params.h_B2 = 0.15;
params.h_P1 = 0.15;
params.h_P2 = 0.15;

% [m] radius of sidearm
params.r_M2 = 0.20;

% startung lengths of the ropes
params.lr10 = 23;
params.lr20 = 13;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimizable parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%

params.M1 = 10000;

% M2 = density*volume = 7860*1.257 approximately 10000
params.M2 = 7860*pi*params.r_M2^2*params.l4;

% I = 0.5*m*r^1
% m = density*volume = 7860*pi*h*r^2
% here thickness h of cable reel = 0.15
params.I_B1 = 0.5*7860*pi*params.h_B1*params.r_B1^4;
params.I_B2 = 0.5*7860*pi*params.h_B2*params.r_B2^4;
params.I_P1 = 0.5*7860*pi*params.h_P1*params.r_P1^4;
params.I_P2 = 0.5*7860*pi*params.h_P2*params.r_P2^4;

params.mu_B1 = 0.2;
params.mu_B2 = 0.2;
params.mu_P1 = 0.2;
params.mu_P2 = 0.2;
