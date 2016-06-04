function [ p ] = pFromParams(params)
% Read optimizable values from params object

p = zeros(10,1);

p(1) = params.M1;
p(2) = params.M2;
p(3) = params.I_B1;
p(4) = params.I_B2;
p(5) = params.I_P1;
p(6) = params.I_P2;
p(7) = params.mu_B1;
p(8) = params.mu_B2;
p(9) = params.mu_P1;
p(10) = params.mu_P2;

end
