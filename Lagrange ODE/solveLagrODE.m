function [sol] = solveLagrODE(params, x_0, t_span, tau_B1, tau_B2)
% Compute the solution of the lagrangean ODE
% Input:
%   params      parameters
%   x_0         initial state of [s,theta,sd,thetad]
%   t_span      timespan for the ode
%   tau_B1      torque function
%   tau_B2      torque function
% Output:
%   sol         solution of ODE

% state x
% x = [s; theta; sd; thetad]
% sd = ds/dt, thetad = d theta/dt
%
% control u
% u = [tau_B1; tau_B2]


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
sol = ode45(@(t,x) rhs(t,x,[tau_B1(t),tau_B2(t)]), t_span, x_0);

end
