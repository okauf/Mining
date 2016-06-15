function [ dTds_dx, dVds_dx, dTdtheta_dx, dVdtheta_dx ] = kinPotEDx( params, s, theta, sd, thetad)
% Calculate derivatives w.r.t. x
% Output:
%   dTds_dx         d(dT/ds)/dx
%   dVds_dx         d(dV/ds)/dx
%   dTdtheta_dx     d(dT/dtheta)/dx
%   dVdtheta_dx     d(dV/dtheta)/dx



[dlr1d2ds, dlr1d2dtheta] = derRopeLength(params, s, theta, sd, thetad);
[dlr1d2ds_dx, dlr1d2dtheta_dx] = derRopeLengthDx(params, s, theta, sd, thetad);


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


% dTds
dTds_dx = [ thetad^2*(params.M1 + 2*params.M2) + dlr1d2ds_dx(1)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2));
    dlr1d2ds_dx(2)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2));
    dlr1d2ds_dx(3)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2));
    dlr1d2ds_dx(4)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2)) + 2*thetad*(2*params.M2*(params.l5 - params.l4/2 + s) + params.M1*(params.l5 + s)) ];

% dVds only depends on theta
dVds_dx = zeros(4,1);
dVds_dx(2) = params.g*cos(theta)*(params.M1 + params.M2);

% dTdtheta
dTdtheta_dx = (params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2))*dlr1d2dtheta_dx;

% dVdtheta
dVdtheta_dx = [ params.g*cos(theta)*(params.M1 + params.M2);
                -params.g*sin(theta)*(params.M2*(params.l5 - params.l4/2 + s) + params.M1*(params.l5 + s));
                0;
                0 ];

end

