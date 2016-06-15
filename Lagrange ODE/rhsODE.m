function [ rhs, rhs_dp, rhs_dx ] = rhsODE( x,u,params )
% Calculate rhs of ODE for given values of x, u and params
% Output:
%   rhs         value of the right hand side
% Optional Output:
%   rhs_dp      derivative of right hand side w.r.t. p
%   rhs_dx      derivative of right hand side w.r.t. x

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

    if nargout > 2
        % The first two components of rhs depend on sd and thetad
        % The third and fourth component are the derivative of A_inv_b
        A_inv_b_dx = lagrddDx(params);
        rhs_dx = [ 0, 0, 1, 0; 0, 0, 0, 1; A_inv_b_dx(x,u)];
    end
end

end
