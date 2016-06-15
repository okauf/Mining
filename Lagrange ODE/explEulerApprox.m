function approxFct = explEulerApprox(params, x_opt, u, mesh)
% Create an approximation function using explicit Euler
% Input:
%   params      fix parameters, not containing optimizable parameters
%   x_opt       optimal reference value, 4xN
%   u           torque values on supporting points of both functions
%   mesh        time steps of discretization
% Output:
%   approxFct   function depending on parameters to calculate approximation

% number of time points
N = length(mesh)+1;

    function [ x, x_dp]  = explEuler(p)
        % Explicit Euler approximation, with derivative if requested

        % p = [M1,M2,I_B1,I_B2,I_P1,I_P2,mu_B1,mu_B2,mu_P1,mu_P2]
        params = pToParams(p,params);

        % x = [s1,      s2,         s3,     ...;
        %     theta1,   theta2,     theta3, ...;
        %     sd1,      sd2,        sd3,    ...;
        %     thetad1,  thetad2,    thetad3,...];
        x = zeros(4,N);
        x(:,1) = x_opt(:,1);    % initial value is the same

        if nargout == 1
            % only calculate approximation x
            for i = 1:(N-1)
                rhs     = rhsODE(x_opt(:,i),u(:,i),params);
                x(:,i+1)  = x_opt(:,i) + mesh(i)*rhs;
            end
        end

        if nargout > 1
            % calculate approximation and its derivative of x w.r.t. p
            % x_{n+1} = x_n + h*rhs
            % dx_{n+1}/dp = h*d(rhs)/dp
            % 
            % dx_1/dp = 0
            x_dp = zeros(4,N,10);

            for i = 1:(N-1)
                [rhs,rhs_dp] = rhsODE(x_opt(:,i),u(:,i),params);
                x(:,i+1)  = x_opt(:,i) + mesh(i)*rhs;
                x_dp(:,i+1,:) = mesh(i)*rhs_dp;
            end

        end
    end

approxFct = @(p) explEuler(p);

end
