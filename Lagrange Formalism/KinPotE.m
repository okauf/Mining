function [ Der_dTdsd, dTds, dVds, Der_dTdthetad, dTdtheta, dVdtheta ] = KinPotE( params, s, theta, sd, thetad, sdd, thetadd )

% two degrees of freedom:
% s denotes distance from attachment point of rope 2 to crossing
% theta denotes angle between side arm and horizontal

syms M1 M2 g I_B1 I_B2 I_P1 I_P2 r_B1 r_B2 r_P1 r_P2 r_M1 r_M2

% M1 = params.m; %mass of load
% M2 = 500; % mass of side arm
% g = 9.81; % gravity constant
% 
% %inertias of side arm and pulleys
% I_B = [];
% I_P1 = [];
% I_P2 = [];
% 
% r_B = []; % radius of basis cable reel
% r_P1 = []; % radius of pulley 1
% r_P2 = []; % radius of pulley 2

[dlr1d2ds, dlr1d2dtheta, Der_dlr1d2dsd, Der_dlr1d2dthetad] = DerRopeLength(params, s, theta, sd, thetad, sdd, thetadd);

% Derivative of T wrt dot(s) and t
Der_dTdsd = (M1 + M2 +  I_B2/r_B2^2 + I_P2/r_P2^2)*sdd ...
    + 0.5*(I_B1/r_B1^2 + I_P1/r_P1^2) * Der_dlr1d2dsd;

%Derivative of T wrt s
dTds = (M1*(s+params.l5) + 2*M2*(s+params.l5 - ...
    0.5*params.l4))* thetad^2 + 0.5*( I_B1/r_B1^2 + I_P1/r_P1^2) * dlr1d2ds;

%Derivative of V wrt s
dVds = (M1+M2)*g*sin(theta);

%Derivative of T wrt dot(theta) and t
Der_dTdthetad = (M1*(s+params.l5)^2 + M2*(0.25*r_M2^2 + 1/12*params.l4^2 + ...
    2*(s+params.l5-params.l4/2)^2))* thetadd ...
    +(2*M1*(s+params.l5) + 4*M2*(s+params.l5 - params.l4/2)) * sd * thetad ... 
    + 0.5*( I_B1/r_B1^2 + I_P1/r_P1^2 ) * Der_dlr1d2dthetad;

%Derivative of T wrt theta
dTdtheta = 0.5*(I_B1/r_B1^2 + I_P1/r_P1^2) * dlr1d2dtheta;

%Derivative of V wrt theta
dVdtheta = (M1 *(s+params.l5) + M2 * (s+params.l5 - params.l4/2))*g*cos(theta);

end

