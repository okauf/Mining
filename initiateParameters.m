% initiation of model parameters


clear params;

% size of excavator

% [m]
params.l1 = 10; % length from A to C
params.l2 = 5; %length from C to D
params.l3 = params.l1 + params.l2;
params.l4 = 10; %length of movable excavator arm
params.l5 = 5; %length from shovel to attachment point of blue pulley

% [kg]
params.m = 10; % mass of load
params.F = params.m *9.81; % weight

params.ang_base = 1/4*pi;

% startung lengths of the ropes
params.lr10_start = 20;
params.lr20_start = 13;

% [N]
params.Fr10 = 0;
params.Fr = @(t) params.Fr10 + t*0.2;

% starting diameter of steel ropes
params.d10 = 0.05;
params.d20 = 0.05;
params.A0 = pi/4*params.d10^2;

params.E_c = 0.2;
params.mu_c = 0.15;

params.r = 0.05;
params.I = 2;
params.F_c = 2;
