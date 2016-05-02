function [ lr1 ] = approxRope1(params, lr1, lr10)
% ROPE1 Approximate length of the stretched rope 1
%   Input:
%       params      Fixed parameters of the excavator model
%       lr1         Previous (stretched) length of rope 1
%       lr10        Current (unstretched) length of rope 1
%   Output:
%       lr1         Length of the stretched rope 1
%
%   This model contains one stretchable rope (rope 1), and a fixed sidearm.
%   When a force is applied to the rope, it gets stretched. This results in a
%   slightly shifted excavator position and results in a slighthly different
%   force pulling on rope 1. This in turn affects the stretching again.
%
%   Instead of calculating the stable position of the excavator with stretched
%   rope 1, this program neglects the change in force and only calculates the
%   stretching by the initial force.

% use initial length of lr1 for approximating the angles
% previous, lr1 is given, so use this instead
%lr1 = params.lr10;

% rope 2 does not stretch in this model
lr2 = params.lr20;

% this angle would depend on the resulting lr1, use the previous lr1 instead
alpha_a = acos(-((lr2-params.l1+params.l5)^2-params.l2^2-(lr1-params.l3)^2)/(2*params.l2*(lr1-params.l3)));

% this angle would depend on lr1, use previous lr1 instead
theta2_a = acos(-((lr1-params.l3)^2-params.l2^2-(lr2-params.l1+params.l5)^2)/(2*params.l2*(lr2-params.l1+params.l5)));

beta_a = pi-alpha_a-theta2_a;

% the angles of the force parallelogramm depend on the state of the excavator
if alpha_a < (pi/2-params.ang_base)
    beta1_a = pi/2 - params.ang_base - alpha_a;
    beta2_a = beta_a - beta1_a;
else
    beta1_a = alpha_a - pi/2 + params.ang_base;
    beta2_a = pi-beta_a-beta1_a;
end

% the force pulling on rope 1
F1 = params.F*sin(beta2_a)/sin(pi-beta1_a-beta2_a);

% the second force pulls on the sidearm, but it is fixed
% F2 = F*sin(beta1_a)/sin(pi-beta1_a-beta2_a);


% solve equation to get new lr1
syms lr1_sol;
eq1 = F1 == params.E_c*(lr1_sol-lr10)/lr10*pi/4*(params.d0*params.mu_c*(lr1_sol-lr10)/lr10 + params.d0)^2;
eq2 = lr1_sol >= 0;

lr1 = double(solve([eq1,eq2],lr1_sol));





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

