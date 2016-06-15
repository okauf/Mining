function [ dTds, dVds, dTdtheta, dVdtheta ] = kinPotE( params, s, theta, sd, thetad)

% two degrees of freedom:
% s denotes distance from attachment point of rope 2 to crossing
% theta denotes angle between side arm and horizontal


[dlr1d2ds, dlr1d2dtheta] = derRopeLength(params, s, theta, sd, thetad);


%Derivative of T wrt s
dTds = (params.M1*(s+params.l5) + 2*params.M2*(s+params.l5 - ...
    0.5*params.l4))* thetad^2 + 0.5*( params.I_B1/params.r_B1^2 + ...
    params.I_P1/params.r_P1^2) * dlr1d2ds;

%Derivative of V wrt s
dVds = (params.M1+params.M2)*params.g*sin(theta);

%Derivative of T wrt theta
dTdtheta = 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2) * dlr1d2dtheta;

%Derivative of V wrt theta
dVdtheta = (params.M1 *(s+params.l5) + params.M2 * (s+params.l5 ...
    - params.l4/2))*params.g*cos(theta);

end

