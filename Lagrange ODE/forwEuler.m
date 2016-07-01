function x = forwEuler(params, x_0, u, mesh)
% Calculate forward Euler approximation of the ODE
% Input:
%   params      fix parameters, not containing optimizable parameters
%   x_0         initial value
%   u           torque values on supporting points of both functions
%   mesh        time steps of discretization
% Output:
%   x           forward Euler approximation of the ODE

% number of time points
N = length(mesh)+1;

% x = [s1,      s2,         s3,     ...;
%     theta1,   theta2,     theta3, ...;
%     sd1,      sd2,        sd3,    ...;
%     thetad1,  thetad2,    thetad3,...];
x = zeros(4,N);
x(:,1) = x_0(:,1);    % initial value is the same

for i = 1:(N-1)
    rhs     = rhsODE(x(:,i),u(:,i),params);
    x(:,i+1)  = x(:,i) + mesh(i)*rhs;
end


end
