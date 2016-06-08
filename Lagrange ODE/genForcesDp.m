function [ Q_s_dp, Q_theta_dp ] = genForcesDp( params, s, theta, sd, thetad, tau_B1, tau_B2)
% Compute derivative of generalized forces w.r.t. parameters
% Input:
%   params      fixed parameters
%   s,theta,sd,thetad
%   tau_B1
%   tau_B2
% Output:
%   Q_s_dp        derivative of generalized force Q_s w.r.t p
%   Q_theta_dp    derivative of generalized force Q_theta w.r.t. p


% for convenience and readability, rename base_ang
phi = params.ang_base;

% d/dt lr1
lr1d = ((s + params.l5)*sd - params.l2*(sd*cos(phi-theta) + (s + params.l5)*sin(phi-theta)*thetad)) / ...
    sqrt(params.l2^2 + (s + params.l5)^2 - 2*params.l2*(s + params.l5)*cos(phi-theta));



% directional unit vector of forces FA1 and FA2
r1 = (s+params.l5)*[cos(theta);sin(theta)] - params.l2*[cos(phi);sin(phi)];
r1 = r1 / norm(r1);
r2 = [cos(theta);sin(theta)];


% rA1 is the position vector of attachement point A1 for force FA1
% drA1ds and drA1dtheta are the derivatives w.r.t. s and theta
drA1ds = [cos(theta);sin(theta)];
drA1dtheta = (s + params.l5)*[-sin(theta);cos(theta)];

% rA2 is the position vector of attachement point A2 for force FA2
% drA2ds and drA2dtheta are the derivatives w.r.t. s and theta
drA2ds = [cos(theta);sin(theta)];
drA2dtheta = s*[-sin(theta);cos(theta)];


% forces FA1 and FA2
%   FA1 = (tau_B1/params.r_B1 - params.mu_B1 * lr1d/params.r_B1 - params.mu_P1 * lr1d/params.r_P1) * r1;
%   FA2 = (tau_B2/params.r_B2 - params.mu_B2 * sd/params.r_B2 - params.mu_P2 * sd/params.r_P2) * r2;
%
% generalized forces
%   Q_s = drA1ds'*FA1 + drA2ds'*FA2;
%   Q_theta = drA1dtheta'*FA1 + drA2dtheta'*FA2;
%
% The generalized forces only depend on the friction coefficients
Q_s_dp = zeros(10,1);
Q_s_dp(7:10) = [ -lr1d/params.r_B1 * drA1ds'*r1;
                 -sd/params.r_B2   * drA2ds'*r2;
                 -lr1d/params.r_P1 * drA1ds'*r1;
                 -sd/params.r_P2   * drA2ds'*r2];

Q_theta_dp = zeros(1,10);
Q_theta_dp(7:10) = [ -lr1d/params.r_B1 * drA1dtheta'*r1;
                     -sd/params.r_B2   * drA2dtheta'*r2;
                     -lr1d/params.r_P1 * drA1dtheta'*r1;
                     -sd/params.r_P2   * drA2dtheta'*r2];


end

