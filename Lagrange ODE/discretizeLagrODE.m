function [sol_ode] = discretizeLagrODE( x, u, params,mesh)
% Approximate solution of the lagrangean ODE using explicit Euler
% state and control combination is feasible if output vector is zero, i.e.
% ode for dynamics of excavator system is (approximately) fulfilled

N = length(mesh);

sol_ode = zeros(4*N,1);

for i = 1:N
    
    % right hand side of ode
    rhs = LagrODE(x(:,i),u(:,i),params);
    
    sol_ode(4*(i-1)+1:4*i) = x(:,i) + mesh(i)*rhs - x(:,i+1);
    
end

end
