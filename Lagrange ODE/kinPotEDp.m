function [ dTds_dp, dVds_dp, dTdtheta_dp, dVdtheta_dp ] = kinPotEDp( params, s, theta, sd, thetad)
% Calculate derivatives w.r.t. p
% Output:
%   dTds_dp         d(dT/ds)/dp
%   dVds_dp         d(dV/ds)/dp
%   dTdtheta_dp     d(dT/dtheta)/dp
%   dVdtheta_dp     d(dV/dtheta)/dp


% rope lengths does not depend on p
[dlr1d2ds, dlr1d2dtheta] = derRopeLength(params, s, theta, sd, thetad);


% Derivative of T wrt s
%   dTds = (params.M1*(s+params.l5) + 2*params.M2*(s+params.l5 - ...
%       0.5*params.l4))* thetad^2 + 0.5*( params.I_B1/params.r_B1^2 + ...
%       params.I_P1/params.r_P1^2) * dlr1d2ds;
%
% Derivative of V wrt s
%   dVds = (params.M1+params.M2)*params.g*sin(theta);
%
% Derivative of T wrt theta
%   dTdtheta = 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2) * dlr1d2dtheta;
%
% Derivative of V wrt theta
%   dVdtheta = (params.M1 *(s+params.l5) + params.M2 * (s+params.l5 ...
%       - params.l4/2))*params.g*cos(theta);


dTds_dp = zeros(10,1);
dTds_dp(1:5) = [ thetad^2*(s + params.l5);
                 2*thetad^2*(s + params.l5 - 0.5*params.l4);
                 dlr1d2ds/(2*params.r_B1^2);
                 0;
                 dlr1d2ds/(2*params.r_P1)];

% dVds only depends on M1 and M2
dVds_dp = zeros(10,1);
dVds_dp(1:2) = [ params.g*sin(theta);
                 params.g*sin(theta)];

% dTdtheta only depends on I_B1 and I_B2
dTdtheta_dp = zeros(10,1);
dTdtheta_dp(3) = dlr1d2dtheta/(2*params.r_B1^2);
dTdtheta_dp(5) = dlr1d2dtheta/(2*params.r_P1^2);

% dVdtheta only depends on M1 and M2
dVdtheta_dp = zeros(10,1);
dVdtheta_dp(1:2) = [ params.g*cos(theta)*(s + params.l5);
                     params.g*cos(theta)*(s + params.l5 - 0.5*params.l4)];

end

