function [ Der_dlr1d2dthetad ] = help2( )

syms l2 l5 s(t) theta(t) phi sd thetad

l_r1 = sqrt(l2^2 + (s(t)+l5)^2 - 2*l2*(s(t)+l5) * cos(phi - theta(t)));

l_r1d2 = diff(l_r1,t)^2;

l_r1d2 = subs(l_r1d2,diff(theta(t),t),thetad);
l_r1d2 = subs(l_r1d2,diff(s(t),t),sd);

Der_dlr1d2dthetad = diff(l_r1d2,thetad);

end

