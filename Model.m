function [ lr1 ] = Model(params)

lr2 = params.lr20;

alpha = @(lr1) acos(-((lr2-params.l1+params.l5)^2-params.l2^2-(lr1-params.l3)^2)/(2*params.l2*(lr1-params.l3)));

theta2 = @(lr1) acos(-((lr1-params.l3)^2-params.l2^2-(lr2-params.l1+params.l5)^2)/(2*params.l2*(lr2-params.l1+params.l5)));

beta = pi-alpha-theta2;

if alpha < (pi/2-params.ang_base)
    beta1 = pi/2 - params.ang_base - alpha;
    beta2 = beta - beta1;
else
    beta1 = alpha - pi/2 + params.ang_base;
    beta2 = pi-beta-beta1;
end

F1 = F*sin(beta2)/sin(pi-beta1-beta2);
% F2 = F*sin(beta1)/sin(pi-beta1-beta2);

A1 = pi/4*d1^2;
sigma1 = F1/A1;

eps1 = (lr1-lr10)/lr10;

E = sigma1/eps1;

eq1 = (d1-d10)/d10;
mu = eq1/eps1;

end

