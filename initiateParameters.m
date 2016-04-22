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
params.m = 1000; % mass of load
params.F = params.m*9.81;

params.ang_base = 1/4*pi;

% startung lengths of the ropes
params.lr10 = 23;
params.lr20 = 13;

% [N]
params.Fr10 = 0;
params.Fr = @(t) params.Fr10 + t*0.2;

% starting diameter [m] and area [m²] of steel ropes
params.d0 = 0.05;
params.A0 = pi/4*params.d0^2;

% modulus of elasticity of steel [N/m²]
params.E_c = 210*10^9;
