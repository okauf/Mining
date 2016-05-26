function [sol] = solveLagrODE(params)
% Compute the solution of the lagrangean ODE

% state x
% x = [s; theta; sd; thetad]
% sd = ds/dt, thetad = d theta/dt
%
% control u
% u = [tau_P1; tau_P2]

% initial value of [s,theta,sd,thetad]
x_0 = [3; -pi/4; 0; 0];

% timespan for the ode
t_span = [0,1];

% torques
k = 1500;   % fixed, such that the system barely moves
%tau_P1 = @(t) -k;
%tau_P2 = @(t) -k;
%tau_P1 = @(t) -k*([t<0.2*pi]*sin(t) - [t>=0.5 && t<0.8]*(t-0.5) + [t>=0.8]*1);
%tau_P2 = @(t) -k*([t<0.2*pi]*sin(t) - [t>=0.5 && t<0.8]*(t-0.5) + [t>=0.8]*1);
tau_P1 = @(t) -k - 30*k*t;
tau_P2 = @(t) -k - k*t;

% From Langrangean formalism:
% A * [sdd; thetadd] = b
% [sdd; thetadd] = A_inv*b
A_inv_b = lagrdd(params);

% right hand side of ode
% d/dt s        =   sd
% d/dt theta    =   thetad
% d/dt [sd;thetad]  =   A_inv*b
rhs = @(t,x,u) [x(3);
                x(4);
                A_inv_b(x,u)];

% solve ode
sol = ode45(@(t,x) rhs(t,x,[tau_P1(t),tau_P2(t)]), t_span, x_0);

end
