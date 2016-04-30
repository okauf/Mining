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

% propertis of the cable reels
params.r_cr = 0.20;                          % [m] radius of cable reel
params.h_cr = 0.15;                          % [m] height of cable reel
params.rho_cr = 7860;                        % [kg/m³] density of steel 
params.V_cr = pi*params.r_cr^2*params.h_cr;  % [m³] volume of cable reel
params.m_cr = params.rho_cr*params.V_cr;     % [kg] mass of cable reel
params.I_cr = 1/12*params.m_cr*(3*params.r_cr^2+params.h_cr^2);   % [kgm²]inertia of cable ree


% modulus of elasticity of steel [N/m²]
params.E_c = 210*10^9;

% poisson's ratio for steel [1]
params.mu_c = 0.28;
