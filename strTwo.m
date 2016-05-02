function [lr1, lr2, sol] = strTwo(params, tau_p1, tau_p2, tspan)
% Compute stretched rope length for movement with two pulleys.
%   Input:
%       params      fix parameters of the model
%       tau_p1      torque of pulley 1
%       tau_p2      torque of pulley 2
%       tspan       time span of movement
%   Output:
%       lr1         stretched length of rope 1
%       lr2         stretched length of rope 2
%       sol         solution from ode
%   The model does not yet contain the effects of the sidearm, thus the movement
%   of the load should only happen on the right side (-> model picture).

% ODE for the unstretched rope lengths depending on the torques
% x represents [lr1,omega1,lr2,omega2] = [lr1,ω₁,lr2,ω₂]
% u represents [tau_p1,tau_p2] = [τ₁,τ₂]
x_0 = [params.lr10 0 params.lr20 0];        % start values of state variables

% ODE:
%   lr1'  = r*ω₁
%   ω₁'   = τ₁/I
%   lr2'  = r*ω₂
%   ω₂    = τ₂/I
rhs = @(t,x,u) [params.r_cr*x(2);
                u(1)/params.I_cr;
                params.r_cr*x(4);
                u(2)/params.I_cr];

sol = ode45( @(t,x) rhs(t,x,[tau_p1(t),tau_p2(t)]), tspan, x_0);

%{
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
%}

% The unstretched lenght is calculated
% Now approximate the stretched length

lr10 = sol.y(1,:);
lr20 = sol.y(3,:);
lr1 = zeros(size(lr10));
lr2 = zeros(size(lr20));
lr1(1) = lr10(1);       % start with the unstretched length
lr2(1) = lr20(1);

%disp(['timesteps: ',num2str(length(lr1))]);
for i = 1:length(lr1)
    disp(['i: ',num2str(i)]);
    % calculate an approximation for the stretched lenght for every time step
    % for the current state use current unstretched length + last stretching
    if i > 1
        lr1(i) = lr10(i) + lr1(i-1) - lr10(i-1);
        lr2(i) = lr20(i) + lr2(i-1) - lr20(i-1);
    end

    [lr1(i),lr2(i)] = approxBothRopes(params,lr1(i),lr10(i),lr2(i),lr20(i),1e-4);
end



end
