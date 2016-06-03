function [ grads, gradt, Der_gradsd, Der_gradtd] = SymComRope( params )

syms s theta sd thetad

phi = params.ang_base;

% l_r1 = sqrt(params.l2^2 + (s+params.l5)^2 - 2*params.l2*(s+params.l5) * cos(phi - theta));% + const;

l_r1d = 0.5*(2*(s+params.l5)*sd - 2*params.l2 * ( sd * cos(phi - theta) + (s+params.l5) * sin(phi - theta) * ...
    thetad))/(sqrt(params.l2^2 + (s+params.l5)^2 - 2*params.l2*(s+params.l5) * cos(phi - theta)));


grads = diff(l_r1d^2,s);

gradt = diff(l_r1d^2,theta);

% gradsd = diff(l_r1d^2,sd);

gradsd = (2*(sd + 5*thetad*sin(theta - pi/4))*((sd*(2*s + 10))/2 - 5*sd*cos(theta - pi/4) + 5*thetad*sin(theta - pi/4)*(s + 5)))/((s + 5)^2 - cos(theta - pi/4)*(10*s + 50) + 25) - ((2*s - 10*cos(theta - pi/4) + 10)*((sd*(2*s + 10))/2 - 5*sd*cos(theta - pi/4) + 5*thetad*sin(theta - pi/4)*(s + 5))^2)/((s + 5)^2 - cos(theta - pi/4)*(10*s + 50) + 25)^2;
 
% manually determined: Derivative of gradsd wrt t
Der_gradsd = [];

% gradtd = diff(l_r1d^2,thetad);

gradtd = (10*sin(theta - pi/4)*(s + 5)*((sd*(2*s + 10))/2 - 5*sd*cos(theta - pi/4) + 5*thetad*sin(theta - pi/4)*(s + 5)))/((s + 5)^2 - cos(theta - pi/4)*(10*s + 50) + 25);
 
% manually determined: Derivative of gradtd wrt t
Der_gradtd = [];

end

