function [ Der_dlr1d2dthetad ] = help2( )

% symbolic computation of d/dt (d/d thetad (d/dt lr1)^2)

syms l2 l5 s(t) theta(t) phi sd thetad_h sd(t) thetad(t) sdd(t) thetadd(t)

l_r1 = sqrt(l2^2 + (s(t)+l5)^2 - 2*l2*(s(t)+l5) * cos(phi - theta(t)));

l_r1d2 = diff(l_r1,t);

l_r1d2 = l_r1d2^2;

l_r1d2 = subs(l_r1d2,diff(theta(t),t),thetad_h);

dlr1d2dthetad = diff(l_r1d2,thetad_h);

dlr1d2dthetad = subs(dlr1d2dthetad,diff(s(t),t),sd(t));
dlr1d2dthetad = subs(dlr1d2dthetad,thetad_h,thetad(t));

Der_dlr1d2dthetad = diff(dlr1d2dthetad,t);

Der_dlr1d2dthetad = subs(Der_dlr1d2dthetad,diff(theta(t),t),thetad);
Der_dlr1d2dthetad = subs(Der_dlr1d2dthetad,diff(s(t),t),sd);
Der_dlr1d2dthetad = subs(Der_dlr1d2dthetad,diff(sd(t),t),sdd);
Der_dlr1d2dthetad = subs(Der_dlr1d2dthetad,diff(thetad(t),t),thetadd);

end

