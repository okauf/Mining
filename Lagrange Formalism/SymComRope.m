function [ grads, gradt, Der_gradsd, Der_gradtd] = SymComRope( params )

syms t s(t) sd(t) sdd(t) theta(t) thetad(t) thetadd(t)
syms t s_h sd_h theta_h thetad_h

phi = params.ang_base;

l_r1 = sqrt(params.l2^2 + (s(t)+params.l5)^2 - 2*params.l2*(s(t)+params.l5) * cos(phi - theta(t)));% + const;


l_r1d = diff(l_r1,t);
% l_r1d2 = 0.5*(2*(s+params.l5)*sd - 2*params.l2 * ( sd * cos(phi - theta) + (s+params.l5) * sin(phi - theta) * ...
%     thetad))/(sqrt(params.l2^2 + (s+params.l5)^2 - 2*params.l2*(s+params.l5) * cos(phi - theta)));

% l_r1d = subs(l_r1d,diff(s(t),t),sd);
% l_r1d = subs(l_r1d,diff(theta(t),t),thetad);

% diff2 = simplify(l_r1d-l_r1d2);

fct = subs(l_r1d^2,diff(s(t),t),sd_h);
fct = subs(fct,diff(theta(t),t),thetad_h);
fct = subs(fct,s(t),s_h);
fct = subs(fct,theta(t),theta_h);

grads = diff(fct,s_h);

gradt = diff(fct,theta_h);

gradsd = diff(fct,sd_h);

gradsd = subs(gradsd,sd_h,sd(t));
gradsd = subs(gradsd,thetad_h,thetad(t));
gradsd = subs(gradsd,s_h,s(t));
gradsd = subs(gradsd,theta_h,theta(t));

Der_gradsd = diff(gradsd,t);

gradtd = diff(fct,thetad_h);

gradtd = subs(gradtd,sd_h,sd(t));
gradtd = subs(gradtd,thetad_h,thetad(t));
gradtd = subs(gradtd,s_h,s(t));
gradtd = subs(gradtd,theta_h,theta(t));

Der_gradtd = diff(gradtd,t);

grads = subs(grads,sd_h,sd);
grads = subs(grads,s_h,s);
grads = subs(grads,theta_h,theta);
grads = subs(grads,thetad_h,thetad);

gradt = subs(gradt,sd_h,sd);
gradt = subs(gradt,s_h,s);
gradt = subs(gradt,theta_h,theta);
gradt = subs(gradt,thetad_h,thetad);

Der_gradsd = subs(Der_gradsd,diff(theta(t),t),thetad);
Der_gradsd = subs(Der_gradsd,diff(s(t),t),sd);
Der_gradsd = subs(Der_gradsd,diff(sd(t),t),sdd);
Der_gradsd = subs(Der_gradsd,diff(thetad(t),t),thetadd);

Der_gradtd = subs(Der_gradtd,diff(theta(t),t),thetad);
Der_gradtd = subs(Der_gradtd,diff(s(t),t),sd);
Der_gradtd = subs(Der_gradtd,diff(sd(t),t),sdd);
Der_gradtd = subs(Der_gradtd,diff(thetad(t),t),thetadd);

end

