function approxFct = rk3Approx(params, x_opt, u, mesh)
% Create an approximation function using Runge-Kutta of order 3
% Input:
%   params      fix parameters, not containing optimizable parameters
%   x_opt       optimal reference value, 4xN
%   u           torque values on supporting points of both functions
%   mesh        time steps of discretization
% Output:
%   approxFct   function depending on parameters to calculate approximation
%
% Butcher-Tableau for RK of order 3 / Simpson's Rule:
% 0   |
% 1/2 | 1/2
%   1 |  -1    2
% -------------------
%     | 1/6  4/6  1/6
%
% The RK method is then explicitely of the form
% k1 = f(x_n)
% k2 = f(x_n + 0.5*h*k1)
% k3 = f(x_n + h*(-k1 +2*k2))
% x_{n+1} = x_n + h*(1/6*k1 + 4/6*k2 + 1/6*k3)

% number of time points
N = length(mesh)+1;

    function [ x, x_dp ] = rungeKutta3(p)
        % Runge Kutta approximation, with derivative if requested

        % p = [M1,M2,I_B1,I_B2,I_P1,I_P2,mu_B1,mu_B2,mu_P1,mu_P2]
        params = pToParams(p,params);

        % x = [s1,      s2,         s3,     ...;
        %     theta1,   theta2,     theta3, ...;
        %     sd1,      sd2,        sd3,    ...;
        %     thetad1,  thetad2,    thetad3,...];
        x = zeros(4,N);
        x(:,1) = x_opt(:,1);    % initial value is the same

        if nargout == 1
            for i = 1:(N-1)
                % Use Runge Kutta of order 3
                % control u is constant within a mesh step
                k1  = rhsODE(x_opt(:,i),u(:,i),params);
                k2  = rhsODE(x_opt(:,i) + 0.5*mesh(i)*k1,u(:,i),params);
                k3  = rhsODE(x_opt(:,i) + mesh(i)*(-k1 + 2*k2),u(:,i),params);
                x(:,i+1)  = x_opt(:,i) + mesh(i)*(1/6*k1 + 4/6*k2 + 1/6*k3);
            end
        end

        if nargout > 1
            % calculate approximation and its derivative w.r.t. p
            % here f(x) = f(x,u,p), but u and p uninteresting for formulas
            % k1 = f(x_n)
            % k2 = f(x_n + 0.5*h*k1)
            % k3 = f(x_n + h*(-k1 +2*k2))
            % x_{n+1} = x_n + h*(1/6*k1 + 4/6*k2 + 1/6*k3)
            %
            % dk1/dp = df(x_n)/dp
            % dk2/dp = df(x_n+0.5*h*k1)/dp + df(...)/dx*0.5*h*dk1/dp
            % dk3/dp = df(x_n+h*(-k1+2*k2))/dp
            %           + df(...)/dx*h*(-dk1/dp + 2*dk2/dp)
            % dx_{n+1}/dp = h*(1/6*dk1/dp + 4/6*dk2/dp + 1/6*dk3/dp)
            x_dp = zeros(4,N,10);

            for i = 1:(N-1)
                % Use Runge Kutta of order 3
                % control u is constant within a mesh step
                [k1,k1_dp]  = rhsODE(x_opt(:,i),u(:,i),params);
                [k2,k2_dp,k2_dx]  = rhsODE(x_opt(:,i) + 0.5*mesh(i)*k1,u(:,i),params);
                [k3,k3_dp,k3_dx]  = rhsODE(x_opt(:,i) + mesh(i)*(-k1 + 2*k2),u(:,i),params);
                % k2_dp and k3_dp do not yet include terms from chain rule
                k2_dp = k2_dp + k2_dx*0.5*mesh(i)*k1_dp; 
                k3_dp = k3_dp + k3_dx*mesh(i)*(-k1_dp + 2*k2_dp);

                x(:,i+1)  = x_opt(:,i) + mesh(i)*(1/6*k1 + 4/6*k2 + 1/6*k3);
                x_dp(:,i+1,:) = mesh(i)*(1/6*k1_dp + 4/6*k2_dp + 1/6*k3_dp);
            end
        end
    end

approxFct = @(p) rungeKutta3(p);

end
