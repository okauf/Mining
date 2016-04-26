function [sol] = unstrOne()
r_cr = 0.20;   % radius of cable reel
h_cr = 0.15;    % height of cable reel
rho_cr = 7860;  % [kg/m³] density of steel for cable reel
V_cr = pi*r_cr^2*h_cr;  % volume of cable reel
m_cr = rho_cr*V_cr;     % mass of cable reel
I_cr = 1/12*m_cr*(3*r_cr^2+h_cr^2)       % inertia of cable reel

% torque of pulley 1
tau_p1  = @(t) I_cr*sin(t);
rhs = @(t,x,u) [r_cr*x(2); u/I_cr]; % ODE for cable 1
t_int = [0 10];      % time interval
x_0 = [23 0];       % start values of state variables

sol = ode45( @(t,x) rhs(t,x,tau_p1(t)), t_int, x_0);



end
