function [ rhs ] = rhsODE( x,u,params )
% Calculate rhs of ODE for given values of x, u and params

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

rhs = [x(3); x(4); A_inv_b(x,u)];

end
