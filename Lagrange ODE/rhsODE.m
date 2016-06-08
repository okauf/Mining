function [ rhs, rhs_dp ] = rhsODE( x,u,params )
% Calculate rhs of ODE for given values of x, u and params
% Output:
%   rhs         value of the right hand side
% Optional Output:
%   rhs_dp      derivative of right hand side w.r.t. p

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


% calculate derivative of rhs if requested
if nargout > 1
    % The first two components of rhs are constant w.r.t. p
    % The third and fourth component are the derivative of A_inv_b
    A_inv_b_dp = lagrddDp(params);
    rhs_dp = [zeros(1,10); zeros(1,10); A_inv_b_dp(x,u)];
end

end
