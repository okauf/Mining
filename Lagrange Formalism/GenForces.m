function [ Q ] = GenForces( params, s, theta, sd, thetad, tau_B1, tau_B2 lr1d )
% Compute generalized forces
% Input:
%   params      fixed parameters
%   s,theta,sd,thetad
%   tau_B1
%   tau_B2
%   lr1d        velocity of rope 1 d(lr1)/dt
% Output:
%   Q           generalized forces Q_s and Q_theta



%<<< some entries in params are not yet defined >>>
% params.r_B1
% params.r_B2
% params.mu_B
% params.mu_P1


% for convenience and readability, rename base_ang
phi = params.base_ang;

% directional unit vector of force FA1
r1 = (s+params.l5)*[cos(theta);sin(theta)] - params.l2*[cos(phi);sin(phi)];
r1 = r1 / norm(r1);

% forces FA1 and FA2
FA1 = (tau_B1/params.r_B1 - params.mu_B * lr1d/params.r_B1 - params.mu_P1 * lr1d/params.r_P1) * r1;

FA2 = (tau_B2/params.r_B2 - params.mu_B * sd/params.r_B2 - params.mu_P2 * sd/params.r_P2)*[cos(theta);sin(theta)];

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

Q = [Q_s;Q_theta];

end

