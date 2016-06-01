function [ grad1, grad2 ] = SetUpODE()

syms lr1 lr2 F l1 l2 l3 l4 l5 ang_base lr10 lr20 E_c mu_c d10 d20

alpha1 = @(lr1,lr2) acos(-((lr2-l1+l5)^2-l2^2-(lr1-l3)^2)/(2*l2*(lr1-l3)));

theta2 = @(lr1,lr2) acos(-((lr1-l3)^2-l2^2-(lr2-l1+l5)^2)/(2*l2*(lr2-l1+l5)));

beta_sum = @(lr1,lr2) pi-alpha1(lr1,lr2)-theta2(lr1,lr2);

% pay attention!
% if alpha < (pi/2-ang_base)
    beta1 = @(lr1,lr2) pi/2 - ang_base - alpha1(lr1,lr2);
    beta2 = @(lr1,lr2) beta_sum(lr1,lr2) - beta1(lr1,lr2);
% else
%     eq4 = beta1 == alpha - pi/2 + ang_base;
%     eq5 = beta2 == pi-beta-beta1;
% end

F1 = @(lr1,lr2,F) F*sin(beta2(lr1,lr2))/sin(pi-beta1(lr1,lr2)-beta2(lr1,lr2));
% F2 = F*sin(beta1)/sin(pi-beta1-beta2);#

eps1 = @(lr1) (lr1-lr10)/lr10;

epsq1 = @(lr1) mu_c * eps1(lr1);

d1 = @(lr1) d10 * (epsq1(lr1) - 1);

A1 = @(lr1) pi/4*d1(lr1)^2;
sigma1 = @(lr1) E_c * eps1(lr1);

F1_2 = @(lr1) sigma1(lr1) * A1(lr1);

grad1 = gradient(F1,[lr1,lr2]);

grad2 = gradient(F1_2,[lr1]);


end

