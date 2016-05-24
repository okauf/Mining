function [ Q_s,Q_theta ] = genForces( params, s, theta, sd, thetad, tau_B1, tau_B2)
% Compute generalized forces
% Input:
%   params      fixed parameters
%   s,theta,sd,thetad
%   tau_B1
%   tau_B2
% Output:
%   Q_s           generalized forces Q_s and Q_theta
%   Q_theta       generalized forces Q_s and Q_theta


% for convenience and readability, rename base_ang
phi = params.ang_base;

% d/dt lr1
lr1d = ((s + params.l5)*sd - params.l2*(sd*cos(phi-theta) + (s + params.l5)*sin(phi-theta)*thetad)) / ...
    sqrt(params.l2^2 + (s + params.l5)^2 - 2*params.l2*(s + params.l5)*cos(phi-theta));



% directional unit vector of force FA1
r1 = (s+params.l5)*[cos(theta);sin(theta)] - params.l2*[cos(phi);sin(phi)];
r1 = r1 / norm(r1);

% forces FA1 and FA2
FA1 = (tau_B1/params.r_B1 - params.mu_B1 * lr1d/params.r_B1 - params.mu_P1 * lr1d/params.r_P1) * r1;

FA2 = (tau_B2/params.r_B2 - params.mu_B2 * sd/params.r_B2 - params.mu_P2 * sd/params.r_P2)*[cos(theta);sin(theta)];

% rA1 is the position vector of attachement point A1 for force FA1
% drA1ds and drA1dtheta are the derivatives w.r.t. s and theta
drA1ds = [cos(theta);sin(theta)];
drA1dtheta = (s + params.l5)*[-sin(theta);cos(theta)];

% rA2 is the position vector of attachement point A2 for force FA2
% drA2ds and drA2dtheta are the derivatives w.r.t. s and theta
drA2ds = [cos(theta);sin(theta)];
drA2dtheta = s*[-sin(theta);cos(theta)];


% generalized forces
Q_s = drA1ds'*FA1 + drA2ds'*FA2;
Q_theta = drA1dtheta'*FA1 + drA2dtheta'*FA2;


end

