% Computation of discrete constraint ensuring validity of Lagrange ODE for
% ODE solution computed with ode45 as input for state and control

initiateParameters;

[sol] = solveLagrODE(params);

%%
% grid for discretization
N = length(sol.x)-1;                  % number of supporting points
int = 1;                              % length time interval
mesh = zeros(1,N);
for i=1:N                             % mesh used for discretization
    mesh(i) = sol.x(i+1)-sol.x(i);
end 

%%
% torques
k = 1500;
tau_P1 = @(t) -k - 30*k*t;
tau_P2 = @(t) -k - k*t;

u = zeros(2,N+1);
u(1,:) = tau_P1(sol.x);
u(2,:) = tau_P2(sol.x);

%%
sol_ode = discretizeLagrODE( sol.y , u, params, mesh);

