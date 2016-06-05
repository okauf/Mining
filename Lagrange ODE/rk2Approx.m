function approxFct = rk2Approx(params, x_opt, u, mesh)
% Create an approximation function using Runge-Kutta of order 2
% Input:
%   params      fix parameters, not containing optimizable parameters
%   x_opt       optimal reference value, 4xN
%   u           torque values on supporting points of both functions
%   mesh        time steps of discretization
% Output:
%   approxFct   function depending on parameters to calculate approximation
%
% Butcher-Tableau for RK of order 2:
% 0   |
% 1/2 | 1/2
% ---------------
%     |   0    1
%
% The RK method is then explicitely of the form
% k1 = f(x_n)
% k2 = f(x_n + 0.5*h*k1)
% x_{n+1} = x_n + h*k2

% number of time points
N = length(mesh)+1;

    function x = rungeKutta2(p)
        % p = [M1,M2,I_B1,I_B2,I_P1,I_P2,mu_B1,mu_B2,mu_P1,mu_P2]
        params = pToParams(p,params);

        % x = [s1,      s2,         s3,     ...;
        %     theta1,   theta2,     theta3, ...;
        %     sd1,      sd2,        sd3,    ...;
        %     thetad1,  thetad2,    thetad3,...];
        x = zeros(4,N);
        x(:,1) = x_opt(:,1);    % initial value is the same

        for i = 1:(N-1)
            % Use Runge Kutta of order 2
            % control u is constant within a mesh step
            k1  = rhsODE(x_opt(:,i),u(:,i),params);
            k2  = rhsODE(x_opt(:,i) + 0.5*mesh(i)*k1,u(:,i),params);
            x(:,i+1)  = x_opt(:,i) + mesh(i)*k2;
        end
    end

approxFct = @(p) rungeKutta2(p);

end
