function [ gradTs,gradVs,Der_gradTsd,gradTt,gradVt,Der_gradTtd] = SymComKinPot( params )

% syms s theta sd thetad
syms M1 M2 g I_B1 I_B2 I_P1 I_P2 r_B1 r_B2 r_P1 r_P2 r_M1 r_M2

% syms t s(t) sd sd(t) sdd theta(t) thetad thetad(t) thetadd

syms t s(t) sd(t) sdd(t) theta(t) thetad(t) thetadd(t)
syms t s_h sd_h theta_h thetad_h

phi = params.ang_base;

%Effects of the load
v1 = sd_h^2 + (s_h+params.l5)^2 * thetad_h^2;
E_kinM1 = 0.5 * M1 * v1;
V_M1 = M1 * g * (s_h+params.l5) * sin(theta_h);

%Effects of the side arm:
v2 = sd_h^2 + (s_h+params.l5-params.l4/2)^2 * thetad_h^2;
I_M2 = 0.25 * M2 * r_M2^2 + 1/12 * M2 * params.l4^2 + M2 * (s_h+params.l5-params.l4/2)^2;
E_kinM2 = 0.5 * M2 * v2;
E_rotM2 = 0.5 * I_M2 * thetad_h^2;
V_M2 = M2 * g * (s_h+params.l5-params.l4/2) * sin(theta_h);

% Variable to describe speed of pulleys:
% l_r1 = sqrt(params.l2^2 + (s+params.l5)^2 - 2*params.l2*(s+params.l5) * cos(phi - theta));% + const;
l_r1d = 0.5*(2*(s_h+params.l5)*sd_h - 2*params.l2 * ( sd_h * cos(phi - theta_h) + (s_h+params.l5) * sin(phi - theta_h) * ...
    thetad_h))/(sqrt(params.l2^2 + (s_h+params.l5)^2 - 2*params.l2*(s_h+params.l5) * cos(phi - theta_h)));

% Effects of the cable reels in the basis:
E_rotB1 = 0.5 * I_B1 * l_r1d^2/r_B1^2; 
E_rotB2 = 0.5 * I_B2 * sd_h^2/r_B2^2; 

% Effects of the pulleys:
E_rotP1 = 0.5 * I_P1 * l_r1d^2/r_P1^2;
E_rotP2 = 0.5 * I_P2 * sd_h^2/r_P2^2;

T = E_kinM1 + E_kinM2 + E_rotM2 + E_rotB1 + E_rotB2 + E_rotP1 + E_rotP2;
V = V_M1 + V_M2;

gradTs = diff(T,s_h);
gradVs = diff(V,s_h);
gradTsd = diff(T,sd_h);
gradTt = diff(T,theta_h);
gradVt = diff(V,theta_h);
gradTtd = diff(T,thetad_h);

gradTsd = subs(gradTsd,s_h,s(t));
gradTsd = subs(gradTsd,sd_h,sd(t));
gradTsd = subs(gradTsd,theta_h,theta(t));
gradTsd = subs(gradTsd,thetad_h,thetad(t));

Der_gradTsd = diff(gradTsd,t);

Der_gradTsd = subs(Der_gradTsd,diff(theta(t),t),thetad);
Der_gradTsd = subs(Der_gradTsd,diff(s(t),t),sd);
Der_gradTsd = subs(Der_gradTsd,diff(sd(t),t),sdd);
Der_gradTsd = subs(Der_gradTsd,diff(thetad(t),t),thetadd);

gradTtd = subs(gradTtd,s_h,s(t));
gradTtd = subs(gradTtd,sd_h,sd(t));
gradTtd = subs(gradTtd,theta_h,theta(t));
gradTtd = subs(gradTtd,thetad_h,thetad(t));

Der_gradTtd = diff(gradTtd,t);

Der_gradTtd = subs(Der_gradTtd,diff(theta(t),t),thetad);
Der_gradTtd = subs(Der_gradTtd,diff(s(t),t),sd);
Der_gradTtd = subs(Der_gradTtd,diff(sd(t),t),sdd);
Der_gradTtd = subs(Der_gradTtd,diff(thetad(t),t),thetadd);

gradTs = subs(gradTs,s_h,s(t));
gradTs = subs(gradTs,sd_h,sd(t));
gradTs = subs(gradTs,theta_h,theta(t));
gradTs = subs(gradTs,thetad_h,thetad(t));

gradVs = subs(gradVs,s_h,s(t));
gradVs = subs(gradVs,sd_h,sd(t));
gradVs = subs(gradVs,theta_h,theta(t));
gradVs = subs(gradVs,thetad_h,thetad(t));

gradTt = subs(gradTt,s_h,s(t));
gradTt = subs(gradTt,sd_h,sd(t));
gradTt = subs(gradTt,theta_h,theta(t));
gradTt = subs(gradTt,thetad_h,thetad(t));

gradVt = subs(gradVt,s_h,s(t));
gradVt = subs(gradVt,sd_h,sd(t));
gradVt = subs(gradVt,theta_h,theta(t));
gradVt = subs(gradVt,thetad_h,thetad(t));


end

