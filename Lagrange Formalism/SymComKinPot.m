function [ gradTs,gradVs,Der_gradTsd,gradTt,gradVt,Der_gradTtd] = SymComKinPot( params )

syms s theta sd thetad
syms M1 M2 g I_B1 I_B2 I_P1 I_P2 r_B1 r_B2 r_P1 r_P2 r_M1 r_M2

% syms t
% s = sym('s(t)');
% theta = sym('theta(t)');
% sd = sym('sd(t)');
% thetad = sym('thetad(t)');

phi = params.ang_base;

%Effects of the load
v1 = sd^2 + (s+params.l5)^2 * thetad^2;
E_kinM1 = 0.5 * M1 * v1;
V_M1 = M1 * g * (s+params.l5) * sin(theta);

%Effects of the side arm:
v2 = sd^2 + (s+params.l5-params.l4/2)^2 * thetad^2;
I_M2 = 0.25 * M2 * r_M2^2 + 1/12 * M2 * params.l4^2 + M2 * (s+params.l5-params.l4/2)^2;
E_kinM2 = 0.5 * M2 * v2;
E_rotM2 = 0.5 * I_M2 * thetad^2;
V_M2 = M2 * g * (s+params.l5-params.l4/2) * sin(theta);

% Variable to describe speed of pulleys:
% l_r1 = sqrt(params.l2^2 + (s+params.l5)^2 - 2*params.l2*(s+params.l5) * cos(phi - theta));% + const;
l_r1d = 0.5*(2*(s+params.l5)*sd - 2*params.l2 * ( sd * cos(phi - theta) + (s+params.l5) * sin(phi - theta) * ...
    thetad))/(sqrt(params.l2^2 + (s+params.l5)^2 - 2*params.l2*(s+params.l5) * cos(phi - theta)));


% Effects of the cable reels in the basis:
E_rotB1 = 0.5 * I_B1 * l_r1d^2/r_B1^2; 
E_rotB2 = 0.5 * I_B2 * sd^2/r_B2^2; 

% Effects of the pulleys:
E_rotP1 = 0.5 * I_P1 * l_r1d^2/r_P1^2;
E_rotP2 = 0.5 * I_P2 * sd^2/r_P2^2;

T = E_kinM1 + E_kinM2 + E_rotM2 + E_rotB1 + E_rotB2 + E_rotP1 + E_rotP2;
V = V_M1 + V_M2;

gradTs = diff(T,s);
gradVs = diff(V,s);

gradTsd = diff(T,sd);

% manually determined: Derivative of gradTsd wrt t
Der_gradTsd = diff(gradTsd,t);

gradTt = diff(T,theta);
gradVt = diff(V,theta);

gradTtd = diff(T,thetad);

% manually determined: Derivative of gradTtd wrt t
Der_gradTtd = diff(gradTtd,t);

end

