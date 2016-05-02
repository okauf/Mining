function [lr1, sol] = strOne(params)
% Compute stretched rope length for movement with one pulley.
%   Input:
%       params      fix parameters of the model
%   Output:
%       lr1         stretched length
%       sol         solution from ode

%{
initiateParameters;
r_cr = 0.20;   % radius of cable reel
h_cr = 0.15;    % height of cable reel
rho_cr = 7860;  % [kg/m³] density of steel for cable reel
V_cr = pi*r_cr^2*h_cr;  % volume of cable reel
m_cr = rho_cr*V_cr;     % mass of cable reel
I_cr = 1/12*m_cr*(3*r_cr^2+h_cr^2);       % inertia of cable reel
%}

% torque of pulley 1
tau_p1  = @(t) params.I_cr*sin(t);


rhs = @(t,x,u) [params.r_cr*x(2); u/params.I_cr]; % ODE for cable 1
t_int = [0 10];      % time interval
x_0 = [params.lr10 0];       % start values of state variables

sol = ode45( @(t,x) rhs(t,x,tau_p1(t)), t_int, x_0);


% The unstretched lenght is calculated
% Now approximate the stretched length

lr10 = sol.y(1,:);
lr1 = zeros(size(lr10));
lr1(1) = lr10(1);   % start with the unstretched length

for i = 1:length(lr1)
    % calculate an approximation for the stretched lenght for every time step
    % use current unstretched length + last stretching for the current state
    if i > 1
        lr1(i) = lr10(i) + lr1(i-1) - lr10(i-1);
    end
    lr1_tmp = lr1(i);    % memory of the last stretch approximation
    for j = 1:10
        % calculate the new approximation
        lr1(i)= approxRope1(params, lr1(i), lr10(i));

        % is this approx good enough? (-> relative error with last approx)
        if abs(lr1(i) - lr1_tmp)/lr1(i) < 1e-4
            break;
        else
            lr1_tmp = lr1(i);
        end
    end
end



end
