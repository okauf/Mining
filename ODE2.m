% function [ grad ] = ODE2( params ) %params

% l1 = params.l1; l2 = params.l2; l3 = params.l3; l4 = params.l4; l5 = params.l5;
% ang_base = params.ang_base; F = params.F; lr10 = params.lr10; lr20 = params.lr20;
% E_c = params.E_c; mu_c = params.mu_c; d10 = params.d10; d20 = params.d20;

syms lr1 lr2
syms F l1 l2 l3 l4 l5 ang_base lr10_start lr20 E_c mu_c d10 d20
syms I r F_c
syms dlr1 dlr2
syms ddlr1 ddlr2

alpha = acos(-((lr2-l1-l5)^2-l2^2-(lr1-l3)^2)/(2*l2*(lr1-l3)));

% derivative of alpha wrt lr1
% alphad = gradient(alpha,lr1);
alphad = ((l2^2 - (l1 + l5 - lr2)^2 + (l3 - lr1)^2)/(2*l2*(l3 - lr1)^2) - (2*l3 - 2*lr1)/(2*l2*(l3 - lr1)))/(1 - (l2^2 - (l1 + l5 - lr2)^2 + (l3 - lr1)^2)^2/(4*l2^2*(l3 - lr1)^2))^(1/2);

% derivative of alpha wrt t
Der_alpha = alphad * dlr1;

% 2nd derivative of alpha wrt lr1
% alphadd = gradient(alphad,lr1);
alphadd = (1/(l2*(l3 - lr1)) + (l2^2 - (l1 + l5 - lr2)^2 + (l3 - lr1)^2)/(l2*(l3 - lr1)^3) - (2*l3 - 2*lr1)/(l2*(l3 - lr1)^2))/(1 - (l2^2 - (l1 + l5 - lr2)^2 + (l3 - lr1)^2)^2/(4*l2^2*(l3 - lr1)^2))^(1/2) + (((l2^2 - (l1 + l5 - lr2)^2 + (l3 - lr1)^2)/(2*l2*(l3 - lr1)^2) - (2*l3 - 2*lr1)/(2*l2*(l3 - lr1)))*((l2^2 - (l1 + l5 - lr2)^2 + (l3 - lr1)^2)^2/(2*l2^2*(l3 - lr1)^3) - ((2*l3 - 2*lr1)*(l2^2 - (l1 + l5 - lr2)^2 + (l3 - lr1)^2))/(2*l2^2*(l3 - lr1)^2)))/(2*(1 - (l2^2 - (l1 + l5 - lr2)^2 + (l3 - lr1)^2)^2/(4*l2^2*(l3 - lr1)^2))^(3/2));

% 2nd derivative of alpha wrt t
Der2_alpha = alphadd * dlr1^2 + alphad * ddlr1;

lr10 = lr10_start + r * alpha;

eq = I/r * Der2_alpha + F_c + E_c*(lr1-lr10)/lr10*pi/4*(d10*mu_c*(lr1-lr10)/lr10+d10)^2;

simplify(solve(eq,ddlr1))



% end

