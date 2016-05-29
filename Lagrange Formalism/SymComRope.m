function [ grads, gradt, gradsd, gradtd] = SymComRope( params )

syms s theta sd thetad

phi = params.ang_base;

% l_r1 = sqrt(params.l2^2 + (s+params.l5)^2 - 2*params.l2*(s+params.l5) * cos(phi - theta));% + const;

l_r1d = 0.5*(2*(s+params.l5)*sd - 2*params.l2 * ( sd * cos(phi - theta) + (s+params.l5) * sin(phi - theta) * ...
    thetad))/(sqrt(params.l2^2 + (s+params.l5)^2 - 2*params.l2*(s+params.l5) * cos(phi - theta)));


grads = diff(l_r1d^2,s);

gradt = diff(l_r1d^2,theta);

gradsd = diff(l_r1d^2,sd);

gradtd = diff(l_r1d^2,thetad);

end

