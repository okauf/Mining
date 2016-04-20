function [ S ] = Model(params)

syms lr1 alpha theta2 beta beta1 beta2 F F1

lr2 = params.lr20;

eq1 = alpha == acos(-((lr2-params.l1+params.l5)^2-params.l2^2-(lr1-params.l3)^2)/(2*params.l2*(lr1-params.l3)));

eq2 = theta2 == acos(-((lr1-params.l3)^2-params.l2^2-(lr2-params.l1+params.l5)^2)/(2*params.l2*(lr2-params.l1+params.l5)));

eq3 = beta == pi-alpha-theta2;

% pay attention!
if alpha < (pi/2-params.ang_base)
    eq4 = beta1 == pi/2 - params.ang_base - alpha;
    eq5 = beta2 == beta - beta1;
else
    eq4 = beta1 == alpha - pi/2 + params.ang_base;
    eq5 = beta2 == pi-beta-beta1;
end

eq6 = F1 == F*sin(beta2)/sin(pi-beta1-beta2);
% F2 = F*sin(beta1)/sin(pi-beta1-beta2);

eq7 = A1 == pi/4*d1^2;
eq8 = sigma1 == F1/A1;

eq9 = eps1 == (lr1-lr10)/lr10;

eq10 = E == sigma1/eps1;

eq11 = epsq1 == (d1-d10)/d10;
eq12 = mu == epsq1/eps1;

S = solve(eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12);

end

