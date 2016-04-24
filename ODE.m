function [ grad2 ] = ODE( params ) %params

l1 = params.l1; l2 = params.l2; l3 = params.l3; l4 = params.l4; l5 = params.l5;
ang_base = params.ang_base; F = params.F; lr10 = params.lr10; lr20 = params.lr20;
E_c = params.E_c; mu_c = params.mu_c; d10 = params.d10; d20 = params.d20;

syms lr1 lr2
% syms F l1 l2 l3 l4 l5 ang_base lr10 lr20 E_c mu_c d10 d20
syms dlr1 dlr2

% dlr1 = 0; dlr2 = 0;

grad = (F*cos(acos(((l5 - l1 + lr2)^2 - (l3 - lr1)^2 + l2^2)/(2*l2*(l5 - l1 + lr2))) - acos(((l3 - lr1)^2 - (l5 - l1 + lr2)^2 + l2^2)/(2*l2*(l3 - lr1))))*sin(acos(((l5 - l1 + lr2)^2 - (l3 - lr1)^2 + l2^2)/(2*l2*(l5 - l1 + lr2))) - ang_base - pi/2)*(((2*dlr2*(l5 - l1 + lr2) + 2*dlr1*(l3 - lr1))/(2*l2*(l5 - l1 + lr2)) - (dlr2*((l5 - l1 + lr2)^2 - (l3 - lr1)^2 + l2^2))/(2*l2*(l5 - l1 + lr2)^2))/(1 - ((l5 - l1 + lr2)^2 - (l3 - lr1)^2 + l2^2)^2/(4*l2^2*(l5 - l1 + lr2)^2))^(1/2) + ((2*dlr2*(l5 - l1 + lr2) + 2*dlr1*(l3 - lr1))/(2*l2*(l3 - lr1)) - (dlr1*((l3 - lr1)^2 - (l5 - l1 + lr2)^2 + l2^2))/(2*l2*(l3 - lr1)^2))/(1 - ((l3 - lr1)^2 - (l5 - l1 + lr2)^2 + l2^2)^2/(4*l2^2*(l3 - lr1)^2))^(1/2)))/sin(acos(((l5 - l1 + lr2)^2 - (l3 - lr1)^2 + l2^2)/(2*l2*(l5 - l1 + lr2))) - acos(((l3 - lr1)^2 - (l5 - l1 + lr2)^2 + l2^2)/(2*l2*(l3 - lr1))))^2 - (F*cos(acos(((l5 - l1 + lr2)^2 - (l3 - lr1)^2 + l2^2)/(2*l2*(l5 - l1 + lr2))) - ang_base - pi/2)*((2*dlr2*(l5 - l1 + lr2) + 2*dlr1*(l3 - lr1))/(2*l2*(l5 - l1 + lr2)) - (dlr2*((l5 - l1 + lr2)^2 - (l3 - lr1)^2 + l2^2))/(2*l2*(l5 - l1 + lr2)^2)))/(sin(acos(((l5 - l1 + lr2)^2 - (l3 - lr1)^2 + l2^2)/(2*l2*(l5 - l1 + lr2))) - acos(((l3 - lr1)^2 - (l5 - l1 + lr2)^2 + l2^2)/(2*l2*(l3 - lr1))))*(1 - ((l5 - l1 + lr2)^2 - (l3 - lr1)^2 + l2^2)^2/(4*l2^2*(l5 - l1 + lr2)^2))^(1/2)) - (E_c*pi*d10^2*dlr1^2*((mu_c*(- lr1*dlr1 + lr10))/lr10 + 1)^2)/(4*lr10) - (E_c*pi*d10^2*dlr1^2*mu_c*((mu_c*(- lr2*dlr1 + lr10))/lr10 + 1)*(- lr1*dlr1 + lr10))/(2*lr10^2);

grad2 = solve(grad == 0,dlr1);

end

