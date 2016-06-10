function [ dlr1d2ds_dx, dlr1d2dtheta_dx] = derRopeLengthDx( params, s, theta, sd, thetad)
% Calculate derivatives from derRopeLength wrt x
% Output:
%   dlr1d2ds_dx     Derivative of dot(lr1)^2 wrt s      wrt x
%   dlr1d2dtheta_dx Derivative of dot(lr1)^2 wrt theta  wrt x

phi = params.ang_base;
l2 = params.l2;
l4 = params.l4;
l5 = params.l5;


dlr1d2ds_dx = zeros(4,1);

dlr1d2ds_dx(1) = (2*((2*l5 + 2*s - 2*l2*cos(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta))^2 + 2*(sd - l2*thetad*sin(phi - theta))*((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)))*(2*l5 + 2*s - 2*l2*cos(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^3 - (2*l2^2*(sd^2*cos(phi - theta)^2 - sd^2 - l2^2*thetad^2*sin(phi - theta)^2 - l5*sd*thetad*sin(2*phi - 2*theta) - s*sd*thetad*sin(2*phi - 2*theta) + 2*l2*l5*thetad^2*(cos(phi - theta) - cos(phi - theta)^3) + 2*l2*s*thetad^2*(cos(phi - theta) - cos(phi - theta)^3) + 2*l2*sd*thetad*sin(phi - theta)))/(2*l5*s + l2^2 + l5^2 + s^2 - 2*l2*l5*cos(phi - theta) - 2*l2*s*cos(phi - theta))^2;



dlr1d2ds_dx(2) = - ((l2^3*sd^2*sin(phi - theta))/2 - (3*l2^3*sd^2*sin(3*phi - 3*theta))/2 - (3*l2^3*l5^2*thetad^2*sin(3*phi - 3*theta))/2 - (3*l2^3*s^2*thetad^2*sin(3*phi - 3*theta))/2 + 2*l2^4*sd*thetad*cos(2*phi - 2*theta) + 2*l2^2*l5*sd^2*sin(2*phi - 2*theta) + 2*l2^4*l5*thetad^2*sin(2*phi - 2*theta) + (l2^3*l5^2*thetad^2*sin(phi - theta))/2 + 2*l2^2*s*sd^2*sin(2*phi - 2*theta) + 2*l2^4*s*thetad^2*sin(2*phi - 2*theta) + (l2^3*s^2*thetad^2*sin(phi - theta))/2 + 2*l2^2*l5^2*sd*thetad*cos(2*phi - 2*theta) + 2*l2^2*s^2*sd*thetad*cos(2*phi - 2*theta) - 3*l2^3*l5*s*thetad^2*sin(3*phi - 3*theta) - 4*l2^3*l5*sd*thetad*cos(phi - theta) - 4*l2^3*s*sd*thetad*cos(phi - theta) + l2^3*l5*s*thetad^2*sin(phi - theta) + 4*l2^2*l5*s*sd*thetad*cos(2*phi - 2*theta))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2 - (4*l2*sin(phi - theta)*(l5 + s)*((2*l5 + 2*s - 2*l2*cos(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta))^2 + 2*(sd - l2*thetad*sin(phi - theta))*((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta))))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^3;



dlr1d2ds_dx(3) = (((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(4*sd*(l5 + s) + l2^2*thetad*sin(2*phi - 2*theta) - 4*l2*sd*cos(phi - theta) - l2*thetad*sin(phi - theta)*(4*l5 + 4*s)) + 2*(l5 + s - l2*cos(phi - theta))*(2*l5 + 2*s - 2*l2*cos(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



dlr1d2ds_dx(4) = -((2*l2^2*thetad*sin(phi - theta)^2*(l5 + s)^2 - 2*l2*sd*sin(phi - theta)*(l5 + s)^2 + 2*l2^2*sd*cos(phi - theta)*sin(phi - theta)*(l5 + s))*(2*l5 + 2*s - 2*l2*cos(phi - theta)) - (4*l2^2*thetad*sin(phi - theta)^2*(l5 + s) - l2*sd*sin(phi - theta)*(4*l5 + 4*s) + 2*l2^2*sd*cos(phi - theta)*sin(phi - theta))*((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;





dlr1d2dtheta_dx = zeros(4,1);

dlr1d2dtheta_dx(1) = - (((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(2*l2*sd^2*sin(phi - theta) + l2^2*thetad^2*sin(2*phi - 2*theta)*(2*l5 + 2*s) + 2*l2^2*sd*thetad*cos(phi - theta)^2 - 2*l2^2*sd*thetad*sin(phi - theta)^2 - 2*l2*sd*thetad*cos(phi - theta)*(2*l5 + 2*s)) - 2*l2*sin(phi - theta)*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta))^2 + 2*l2*(2*l5 + 2*s - 2*l2*cos(phi - theta))*(l5*thetad*cos(phi - theta) - sd*sin(phi - theta) + s*thetad*cos(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)) + 4*l2*sin(phi - theta)*(sd - l2*thetad*sin(phi - theta))*(l5 + s)*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2 - (2*(2*l2*sin(phi - theta)*(l5 + s)*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta))^2 - 2*l2*((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(l5*thetad*cos(phi - theta) - sd*sin(phi - theta) + s*thetad*cos(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)))*(2*l5 + 2*s - 2*l2*cos(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^3;



dlr1d2dtheta_dx(2) = (((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(2*l2^2*sd^2*sin(phi - theta)^2 - 2*l2^2*sd^2*cos(phi - theta)^2 + 2*l2^2*thetad^2*cos(phi - theta)^2*(l5 + s)^2 - 2*l2^2*thetad^2*sin(phi - theta)^2*(l5 + s)^2 + l2*sd^2*cos(phi - theta)*(2*l5 + 2*s) + 2*l2*sd*thetad*sin(phi - theta)*(l5 + s)^2 - 4*l2^2*sd*thetad*sin(2*phi - 2*theta)*(l5 + s)) - 2*l2*cos(phi - theta)*(l5 + s)*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta))^2)/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2 + (4*l2*sin(phi - theta)*(l5 + s)*(2*l2*sin(phi - theta)*(l5 + s)*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta))^2 - 2*l2*((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(l5*thetad*cos(phi - theta) - sd*sin(phi - theta) + s*thetad*cos(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta))))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^3;



dlr1d2dtheta_dx(3) = (((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(2*l2^2*sd*sin(2*phi - 2*theta) - 2*l2^2*thetad*cos(phi - theta)^2*(l5 + s) - 2*l2*sd*sin(phi - theta)*(2*l5 + 2*s) + 2*l2^2*thetad*sin(phi - theta)^2*(l5 + s) + 2*l2*thetad*cos(phi - theta)*(l5 + s)^2) - 4*l2*sin(phi - theta)*(l5 + s)*(l5 + s - l2*cos(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



dlr1d2dtheta_dx(4) = -(((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(2*l2^2*sd*cos(phi - theta)^2*(l5 + s) - 2*l2^2*sd*sin(phi - theta)^2*(l5 + s) - 2*l2*sd*cos(phi - theta)*(l5 + s)^2 + 4*l2^2*thetad*cos(phi - theta)*sin(phi - theta)*(l5 + s)^2) - 2*l2*sin(phi - theta)*(l5 + s)*(2*l2^2*thetad*sin(phi - theta)^2*(l5 + s)^2 - 2*l2*sd*sin(phi - theta)*(l5 + s)^2 + 2*l2^2*sd*cos(phi - theta)*sin(phi - theta)*(l5 + s)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



end
