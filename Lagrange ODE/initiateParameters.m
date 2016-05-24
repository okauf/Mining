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
params.r_B1 = 0.20;
params.r_B2 = 0.20;
params.r_P1 = 0.20;
params.r_P2 = 0.20;

% [m] radius of sidearm
params.r_M2 = 0.10;

% startung lengths of the ropes
params.lr10 = 23;
params.lr20 = 13;

% optimizable parameters
% not yet reasonable values
params.M1 = 1000;
params.M2 = 1000;
params.I_B1 = 100;
params.I_B2 = 100;
params.I_P1 = 100;
params.I_P2 = 100;
params.mu_B1 = 10;
params.mu_B2 = 10;
params.mu_P1 = 10;
params.mu_P2 = 10;
