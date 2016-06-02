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

    function x = explEuler(p)
        % p = [M1,M2,I_B1,I_B2,I_P1,I_P2,mu_B1,mu_B2,mu_P1,mu_P2]
        params.M1 = p(1);
        params.M2 = p(2);
        params.I_B1 = p(3);
        params.I_B2 = p(4);
        params.I_P1 = p(5);
        params.I_P2 = p(6);
        params.mu_B1 = p(7);
        params.mu_B2 = p(8);
        params.mu_P1 = p(9);
        params.mu_P2 = p(10);

        % x = [s1,      s2,         s3,     ...;
        %     theta1,   theta2,     theta3, ...;
        %     sd1,      sd2,        sd3,    ...;
        %     thetad1,  thetad2,    thetad3,...];
        x = zeros(4,N);
        x(:,1) = x_opt(:,1);    % initial value is the same

        for i = 1:(N-1)
            rhs     = rhsODE(x_opt(:,i),u(:,i),params);
            x(:,i+1)  = x_opt(:,i) + mesh(i)*rhs;
        end
    end

approxFct = @(p) explEuler(p);

end
