function [ lr1 ] = approx_rope1(params, lr1)
% ROPE1 Approximate length of the stretched rope 1
%   Input:
%       params      Fixed parameters of the excavator model
%       lr1         Current (stretched) length of lr1
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
% currently, lr1 is given, so use this instead
%lr1 = params.lr10;

% rope 2 does not stretch in this model
lr2 = params.lr20;

% this angle would depend on lr1, use params.lr1 instead
alpha_a = acos(-((lr2-params.l1+params.l5)^2-params.l2^2-(lr1-params.l3)^2)/(2*params.l2*(lr1-params.l3)));

% this angle would depend on lr1, use params.lr1 instead
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

% tensile stress. [N/m²]
% instead of using the area after the stretching, use the initial area
sigma1 = F1/params.A0;

% extensin [1]
eps1 = sigma1/params.E_c;

lr1 = params.lr10*(eps1 + 1);

disp(['alpha_a   ',num2str(alpha_a*180/pi)])
disp(['beta_a    ',num2str(beta_a*180/pi)])
disp(['beta1_a   ',num2str(beta1_a*180/pi)])
disp(['beta2_a   ',num2str(beta2_a*180/pi)])
disp(['theta2_a  ',num2str(theta2_a*180/pi)])
disp(['F1        ',num2str(F1)])
disp(['sigma1    ',num2str(sigma1)])
disp(['eps1      ',num2str(eps1)])
disp(['lr1       ',num2str(lr1)])

%{
% No feedback loop, so following equations were not used as usually
A1 = pi/4*d1^2;
eps1 = (lr1-lr10)/lr10;
E_c = sigma1/eps1;
epsq1 = (d1-d10)/d10;
mu = epsq1/eps1;
%}

end

