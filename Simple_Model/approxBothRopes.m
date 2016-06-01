function [ lr1,lr2 ] = approxBothRopes(params, lr1, lr10, lr2, lr20, err)
% Approximate length of both stretched ropes
%   Input:
%       params      Fixed parameters of the excavator model
%       lr1         Previous (stretched) length of rope 1
%       lr10        Current (unstretched) length of rope 1
%       lr2         Previous (stretched) length of rope 2
%       lr20        Current (unstretched) length of rope 2
%       err         error threshold for recursive call
%   Output:
%       lr1         Length of the stretched rope 1
%       lr2         Length of the stretched rope 2
%
%   This model contains both stretchable ropes.
%   When a force is applied to the rope, it gets stretched. This results in a
%   slightly shifted excavator position and results in a slighthly different
%   force pulling on the ropes. This in turn affects the stretching again.
%
%   This model only considers a force from the load weight. It explicitely
%   neglects rotational forces and forces resulting from the sidearm weight.
%
%   Instead of calculating the stable position of the excavator with stretched
%   ropes, this program neglects the change in force and only calculates the
%   stretching by the initial force. As long as the absolute error of the current
%   solution compared to the input is large (>err), recursively repeat until
%   the error threshold is fullfilled.



% this angle would depend on the resulting lr1, use the previous lr1 instead
alpha_a = acos(-((lr2-params.l1+params.l5)^2-params.l2^2-(lr1-params.l3)^2)/(2*params.l2*(lr1-params.l3)));

% this angle would depend on lr1, use previous lr1 instead
theta2_a = acos(-((lr1-params.l3)^2-params.l2^2-(lr2-params.l1+params.l5)^2)/(2*params.l2*(lr2-params.l1+params.l5)));

beta_a = pi-alpha_a-theta2_a;

% the angles of the force parallelogramm depend on the state of the excavator
if alpha_a < (pi/2-params.ang_base)
    % the load is on the left side
    beta1_a = pi/2 - params.ang_base - alpha_a;
    beta2_a = beta_a - beta1_a;
else
    % the load is on the right side
    beta1_a = alpha_a - pi/2 + params.ang_base;
    beta2_a = pi-beta_a-beta1_a;
    disp('The load is on the right side, but sidearm mass not considered yet.');
end

% the force pulling in direction of rope 1
F1 = params.F*sin(beta2_a)/sin(pi-beta1_a-beta2_a);

% the force pulling in direction of the sidearm / rope 2
F2 = params.F*sin(beta1_a)/sin(pi-beta1_a-beta2_a);


% solve equation to get new lr1
syms lr1_sol;
eq1 = F1 == params.E_c*(lr1_sol-lr10)/lr10*pi/4*(params.d0*params.mu_c*(lr1_sol-lr10)/lr10 + params.d0)^2;
eq2 = lr1_sol >= 0;

lr1_tmp = double(solve([eq1,eq2],lr1_sol));

% solve equation to get new lr2
syms lr2_sol;
eq3 = F2 == params.E_c*(lr2_sol-lr20)/lr20*pi/4*(params.d0*params.mu_c*(lr2_sol-lr20)/lr20 + params.d0)^2;
eq4 = lr2_sol >= 0;

lr2_tmp = double(solve([eq3,eq4],lr2_sol));


% Repeat the approximation recursively, if the absolute error is too large
abs_err1 = abs(lr1 - lr1_tmp);
abs_err2 = abs(lr2 - lr2_tmp);
%disp(['lr1: ',num2str(lr1),', lr2: ',num2str(lr2)]);
%disp(['lr1_tmp: ',num2str(lr1_tmp),', lr2_tmp: ',num2str(lr2_tmp)]);
%disp(['F1: ',num2str(F1),', F2: ',num2str(F2)]);
%disp(['err1: ',num2str(abs_err1), ', err2: ',num2str(abs_err2)]);

if abs_err1 > err || abs_err2 > err
    disp('approxBothRopes: error too large, repeat!');
    % >>>
    if abs(lr1_tmp) > 40
        disp('...complex, returning')
        return;
    end
    % <<<
    [lr1_tmp,lr2_tmp] = approxBothRopes(params,lr1_tmp,lr10,lr2_tmp,lr20,err);
end

lr1 = lr1_tmp;
lr2 = lr2_tmp;


%{
disp(['alpha_a   ',num2str(alpha_a*180/pi)])
disp(['beta_a    ',num2str(beta_a*180/pi)])
disp(['beta1_a   ',num2str(beta1_a*180/pi)])
disp(['beta2_a   ',num2str(beta2_a*180/pi)])
disp(['theta2_a  ',num2str(theta2_a*180/pi)])
disp(['F1        ',num2str(F1)])
disp(['sigma1    ',num2str(sigma1)])
disp(['eps1      ',num2str(eps1)])
disp(['lr1       ',num2str(lr1)])
%}


end

