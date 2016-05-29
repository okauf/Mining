function [ rhs ] = LagrODE( x,u,params ) % further required output: Derivative rhs wrt params (optimization var)

% right hand side of ode

A_inv_b = lagrdd(params);

rhs = [x(3); x(4); A_inv_b(x,u)];  % A\b

end

