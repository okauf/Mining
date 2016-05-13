function [ Q ] = GenForces( params )

phi = params.base_ang;

r1 = (s+params.l5)*[cos(theta);sin(theta)] - params.l2*[cos(phi);sin(phi)];
r1 = r1 / norm(r1);

F1 = (tau_B1/r_B1 - mu_B * lr1d/r_B1 - mu_P1 * lr1d/r_P1) * r1;

F2 = (tau_B2/r_B2 - mu_B * sd/r_B2 - mu_P2 * sd/r_P2)*[cos(theta);sin(theta)];

Q = (F1 + F2); %XXX

end

