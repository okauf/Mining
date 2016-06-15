function [a_s_dx,b_s_dx,c_s_dx,a_theta_dx,b_theta_dx,c_theta_dx] = derderRopeLengthDx(params, s, theta, sd, thetad)
% Calculate derivatives from derderRopeLength wrt x

phi = params.ang_base;
l2 = params.l2;
l4 = params.l4;
l5 = params.l5;


% a_s
a_s_dx = zeros(4,1);

a_s_dx(1) = -(4*l2^2*(cos(phi - theta)^2 - 1)*(l5 + s - l2*cos(phi - theta)))/(2*l5*s + l2^2 + l5^2 + s^2 - 2*l2*l5*cos(phi - theta) - 2*l2*s*cos(phi - theta))^2;



a_s_dx(2) = (4*cos(phi - theta)*sin(phi - theta)*l2^2 - sin(phi - theta)*(4*l5 + 4*s)*l2)/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s)) + (2*l2*sin(phi - theta)*(l5 + s)*(2*l2^2*cos(phi - theta)^2 + 2*(l5 + s)^2 - l2*cos(phi - theta)*(4*l5 + 4*s)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



% b_s
b_s_dx = zeros(4,1);

b_s_dx(1) = (sin(2*phi - 2*theta)*l2^2 - 2*sin(phi - theta)*(2*l5 + 2*s)*l2)/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s)) + ((2*l2*sin(phi - theta)*(l5 + s)^2 - l2^2*sin(2*phi - 2*theta)*(l5 + s))*(2*l5 + 2*s - 2*l2*cos(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



b_s_dx(2) = (2*l2*cos(phi - theta)*(l5 + s)^2 - 2*l2^2*cos(phi - theta)^2*(l5 + s) + 2*l2^2*sin(phi - theta)^2*(l5 + s))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s)) - (2*l2*sin(phi - theta)*(l5 + s)*(2*l2*sin(phi - theta)*(l5 + s)^2 - 2*l2^2*cos(phi - theta)*sin(phi - theta)*(l5 + s)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;




% c_s
c_s_dx = zeros(4,1);

c_s_dx(1) = (2*l2^2*sd^2 - 2*l2^2*sd^2*cos(2*phi - 2*theta) - 2*l2^4*thetad^2*cos(2*phi - 2*theta) - 12*l2^2*l5^2*thetad^2 - 12*l2^2*s^2*thetad^2 - 9*l2^3*sd*thetad*sin(phi - theta) - 6*l2^2*l5^2*thetad^2*cos(2*phi - 2*theta) - 6*l2^2*s^2*thetad^2*cos(2*phi - 2*theta) + 8*l2*l5^3*thetad^2*cos(phi - theta) + 10*l2^3*l5*thetad^2*cos(phi - theta) + 8*l2*s^3*thetad^2*cos(phi - theta) + 10*l2^3*s*thetad^2*cos(phi - theta) - l2^3*sd*thetad*sin(3*phi - 3*theta) + 2*l2^3*l5*thetad^2*cos(3*phi - 3*theta) + 2*l2^3*s*thetad^2*cos(3*phi - 3*theta) - 24*l2^2*l5*s*thetad^2 - 12*l2^2*l5*s*thetad^2*cos(2*phi - 2*theta) + 24*l2*l5*s^2*thetad^2*cos(phi - theta) + 24*l2*l5^2*s*thetad^2*cos(phi - theta) + 6*l2^2*l5*sd*thetad*sin(2*phi - 2*theta) + 6*l2^2*s*sd*thetad*sin(2*phi - 2*theta))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2 - (2*(((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(4*sd^2*(l5 + s) - 4*l2*sd^2*cos(phi - theta) + 2*l2*thetad^2*cos(phi - theta)*(l5 + s)^2 + 3*l2^2*sd*thetad*sin(2*phi - 2*theta) - 2*l2^2*thetad^2*cos(phi - theta)^2*(l5 + s) + 2*l2^2*thetad^2*sin(phi - theta)^2*(l5 + s) - 8*l2*sd*thetad*sin(phi - theta)*(l5 + s)) - 2*(l5 + s - l2*cos(phi - theta))*(2*l2*sd*cos(phi - theta) - 2*sd*(l5 + s) + 2*l2*thetad*sin(phi - theta)*(l5 + s))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)))*(2*l5 + 2*s - 2*l2*cos(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^3;



c_s_dx(2) = ((2*l2*sd*cos(phi - theta) - 2*sd*(l5 + s) + 2*l2*thetad*sin(phi - theta)*(l5 + s))*(2*l2^2*sd*sin(2*phi - 2*theta) - 2*l2^2*thetad*cos(phi - theta)^2*(l5 + s) + 2*l2^2*thetad*sin(phi - theta)^2*(l5 + s) - 4*l2*sd*sin(phi - theta)*(l5 + s) + 2*l2*thetad*cos(phi - theta)*(l5 + s)^2) - ((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(4*l2*sd^2*sin(phi - theta) - 2*l2*thetad^2*sin(phi - theta)*(l5 + s)^2 + 6*l2^2*sd*thetad*cos(phi - theta)^2 - 6*l2^2*sd*thetad*sin(phi - theta)^2 + 4*l2^2*thetad^2*sin(2*phi - 2*theta)*(l5 + s) - 8*l2*sd*thetad*cos(phi - theta)*(l5 + s)) + 4*l2*(l5 + s - l2*cos(phi - theta))*(l5*thetad*cos(phi - theta) - sd*sin(phi - theta) + s*thetad*cos(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)) - 2*l2*sin(phi - theta)*(l5 + s)*(4*sd^2*(l5 + s) - 4*l2*sd^2*cos(phi - theta) + 2*l2*thetad^2*cos(phi - theta)*(l5 + s)^2 + 3*l2^2*sd*thetad*sin(2*phi - 2*theta) - 2*l2^2*thetad^2*cos(phi - theta)^2*(l5 + s) + 2*l2^2*thetad^2*sin(phi - theta)^2*(l5 + s) - 8*l2*sd*thetad*sin(phi - theta)*(l5 + s)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2 + (4*l2*sin(phi - theta)*(l5 + s)*(((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(4*sd^2*(l5 + s) - 4*l2*sd^2*cos(phi - theta) + 2*l2*thetad^2*cos(phi - theta)*(l5 + s)^2 + 3*l2^2*sd*thetad*sin(2*phi - 2*theta) - 2*l2^2*thetad^2*cos(phi - theta)^2*(l5 + s) + 2*l2^2*thetad^2*sin(phi - theta)^2*(l5 + s) - 8*l2*sd*thetad*sin(phi - theta)*(l5 + s)) - 2*(l5 + s - l2*cos(phi - theta))*(2*l2*sd*cos(phi - theta) - 2*sd*(l5 + s) + 2*l2*thetad*sin(phi - theta)*(l5 + s))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta))))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^3;



c_s_dx(3) = (2*(l5 + s - l2*cos(phi - theta))^2*(2*l2*sd*cos(phi - theta) - 2*sd*(l5 + s) + 2*l2*thetad*sin(phi - theta)*(l5 + s)) + ((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(8*sd*(l5 + s) + 3*l2^2*thetad*sin(2*phi - 2*theta) - 8*l2*sd*cos(phi - theta) - 8*l2*thetad*sin(phi - theta)*(l5 + s)) + 2*(l5 + s - l2*cos(phi - theta))*(2*l5 + 2*s - 2*l2*cos(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



c_s_dx(4) = -(l2*(l5 + s)*(2*l2*sd*cos(phi - theta) - 2*sd*(l5 + s) + 2*l2*thetad*sin(phi - theta)*(l5 + s))*(2*l5*sin(phi - theta) - l2*sin(2*phi - 2*theta) + 2*s*sin(phi - theta)) - ((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(3*l2^2*sd*sin(2*phi - 2*theta) - 4*l2^2*thetad*cos(phi - theta)^2*(l5 + s) + 4*l2^2*thetad*sin(phi - theta)^2*(l5 + s) - 8*l2*sd*sin(phi - theta)*(l5 + s) + 4*l2*thetad*cos(phi - theta)*(l5 + s)^2) + 4*l2*sin(phi - theta)*(l5 + s)*(l5 + s - l2*cos(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



% a_theta
a_theta_dx = zeros(4,1);

a_theta_dx(1) = ((2*sin(phi - theta)*(l5 + s)^2 - l2^2*sin(2*phi - 2*theta)*(l5 + s))*(2*l5 + 2*s - 2*l2*cos(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2 - (- sin(2*phi - 2*theta)*l2^2 + sin(phi - theta)*(4*l5 + 4*s))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s));



a_theta_dx(2) = (2*cos(phi - theta)*(l5 + s)^2 - 2*l2^2*cos(phi - theta)^2*(l5 + s) + 2*l2^2*sin(phi - theta)^2*(l5 + s))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s)) - (2*l2*sin(phi - theta)*(l5 + s)*(2*sin(phi - theta)*(l5 + s)^2 - 2*l2^2*cos(phi - theta)*sin(phi - theta)*(l5 + s)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



% b_theta
b_theta_dx = zeros(4,1);

b_theta_dx(1) = (2*l2^2*sin(phi - theta)^2*(2*l5 + 2*s))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s)) - (2*l2^2*sin(phi - theta)^2*(l5 + s)^2*(2*l5 + 2*s - 2*l2*cos(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



b_theta_dx(2) = (4*l2^3*sin(phi - theta)^3*(l5 + s)^3)/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2 - (4*l2^2*cos(phi - theta)*sin(phi - theta)*(l5 + s)^2)/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s));



% c_theta
c_theta_dx = zeros(4,1);

c_theta_dx(1) = ((2*l5 + 2*s - 2*l2*cos(phi - theta))*(l2^2*sd^2*sin(2*phi - 2*theta) - l2^2*thetad^2*sin(2*phi - 2*theta)*(l5 + s)^2 - l2*sd^2*sin(phi - theta)*(4*l5 + 4*s) - 2*l2^2*sd*thetad*cos(phi - theta)^2*(l5 + s) + 6*l2^2*sd*thetad*sin(phi - theta)^2*(l5 + s) + l2*s*sd*thetad*cos(phi - theta)*(l5 + s)^2) - ((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(4*l2*sd^2*sin(phi - theta) + l2^2*thetad^2*sin(2*phi - 2*theta)*(2*l5 + 2*s) + 2*l2^2*sd*thetad*cos(phi - theta)^2 - 6*l2^2*sd*thetad*sin(phi - theta)^2 - l2*sd*thetad*cos(phi - theta)*(l5 + s)^2 - l2*s*sd*thetad*cos(phi - theta)*(2*l5 + 2*s)) + 2*l2*sin(phi - theta)*(2*l2*sd*cos(phi - theta) - sd*(2*l5 + 2*s) + 2*l2*thetad*sin(phi - theta)*(l5 + s))*(l2*sd*cos(phi - theta) - 2*s*sd - 2*l5*sd + 2*l2*l5*thetad*sin(phi - theta) + 2*l2*s*thetad*sin(phi - theta)) - 2*l2*sin(phi - theta)*(l5 + s)*(2*sd - 2*l2*thetad*sin(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2 - (2*(((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(l2^2*sd^2*sin(2*phi - 2*theta) - l2^2*thetad^2*sin(2*phi - 2*theta)*(l5 + s)^2 - l2*sd^2*sin(phi - theta)*(4*l5 + 4*s) - 2*l2^2*sd*thetad*cos(phi - theta)^2*(l5 + s) + 6*l2^2*sd*thetad*sin(phi - theta)^2*(l5 + s) + l2*s*sd*thetad*cos(phi - theta)*(l5 + s)^2) + 2*l2*sin(phi - theta)*(l5 + s)*(2*l2*sd*cos(phi - theta) - sd*(2*l5 + 2*s) + 2*l2*thetad*sin(phi - theta)*(l5 + s))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)))*(2*l5 + 2*s - 2*l2*cos(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^3;



c_theta_dx(2) = (4*l2*sin(phi - theta)*(l5 + s)*(((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(l2^2*sd^2*sin(2*phi - 2*theta) - l2^2*thetad^2*sin(2*phi - 2*theta)*(l5 + s)^2 - l2*sd^2*sin(phi - theta)*(4*l5 + 4*s) - 2*l2^2*sd*thetad*cos(phi - theta)^2*(l5 + s) + 6*l2^2*sd*thetad*sin(phi - theta)^2*(l5 + s) + l2*s*sd*thetad*cos(phi - theta)*(l5 + s)^2) + 2*l2*sin(phi - theta)*(l5 + s)*(2*l2*sd*cos(phi - theta) - sd*(2*l5 + 2*s) + 2*l2*thetad*sin(phi - theta)*(l5 + s))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta))))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^3 - ((2*l2*sd*cos(phi - theta) - sd*(2*l5 + 2*s) + 2*l2*thetad*sin(phi - theta)*(l5 + s))*(2*l2^2*sd*cos(phi - theta)^2*(l5 + s) - 2*l2^2*sd*sin(phi - theta)^2*(l5 + s) + 2*l2^2*thetad*sin(2*phi - 2*theta)*(l5 + s)^2 - 2*l2*sd*cos(phi - theta)*(l5 + s)^2) - ((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(2*l2^2*sd^2*sin(phi - theta)^2 - 2*l2^2*sd^2*cos(phi - theta)^2 + 2*l2^2*thetad^2*cos(phi - theta)^2*(l5 + s)^2 - 2*l2^2*thetad^2*sin(phi - theta)^2*(l5 + s)^2 + l2*sd^2*cos(phi - theta)*(4*l5 + 4*s) - 8*l2^2*sd*thetad*sin(2*phi - 2*theta)*(l5 + s) + l2*s*sd*thetad*sin(phi - theta)*(l5 + s)^2) + 2*l2*sin(phi - theta)*(l5 + s)*(l2^2*sd^2*sin(2*phi - 2*theta) - l2^2*thetad^2*sin(2*phi - 2*theta)*(l5 + s)^2 - l2*sd^2*sin(phi - theta)*(4*l5 + 4*s) - 2*l2^2*sd*thetad*cos(phi - theta)^2*(l5 + s) + 6*l2^2*sd*thetad*sin(phi - theta)^2*(l5 + s) + l2*s*sd*thetad*cos(phi - theta)*(l5 + s)^2) + 4*l2^2*sin(phi - theta)*(l5 + s)*(l5*thetad*cos(phi - theta) - sd*sin(phi - theta) + s*thetad*cos(phi - theta))*(l2*sd*cos(phi - theta) - s*sd - l5*sd + l2*l5*thetad*sin(phi - theta) + l2*s*thetad*sin(phi - theta)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



c_theta_dx(3) = -((2*l2^2*thetad*sin(phi - theta)^2*(l5 + s)^2 - 2*l2*sd*sin(phi - theta)*(l5 + s)^2 + 2*l2^2*sd*cos(phi - theta)*sin(phi - theta)*(l5 + s))*(2*l5 + 2*s - 2*l2*cos(phi - theta)) - ((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(6*l2^2*thetad*sin(phi - theta)^2*(l5 + s) - 2*l2*sd*sin(phi - theta)*(4*l5 + 4*s) - 2*l2^2*thetad*cos(phi - theta)^2*(l5 + s) + 4*l2^2*sd*cos(phi - theta)*sin(phi - theta) + l2*s*thetad*cos(phi - theta)*(l5 + s)^2) + (2*l2*sin(phi - theta)*(l5 + s)^2 - 2*l2^2*cos(phi - theta)*sin(phi - theta)*(l5 + s))*(2*l2*sd*cos(phi - theta) - sd*(2*l5 + 2*s) + 2*l2*thetad*sin(phi - theta)*(l5 + s)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



c_theta_dx(4) = (2*l2*sin(phi - theta)*(l5 + s)*(2*l2^2*thetad*sin(phi - theta)^2*(l5 + s)^2 - 2*l2*sd*sin(phi - theta)*(l5 + s)^2 + 2*l2^2*sd*cos(phi - theta)*sin(phi - theta)*(l5 + s)) - ((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))*(2*l2^2*sd*cos(phi - theta)^2*(l5 + s) - 6*l2^2*sd*sin(phi - theta)^2*(l5 + s) - l2*s*sd*cos(phi - theta)*(l5 + s)^2 + 4*l2^2*thetad*cos(phi - theta)*sin(phi - theta)*(l5 + s)^2) + 2*l2^2*sin(phi - theta)^2*(l5 + s)^2*(2*l2*sd*cos(phi - theta) - sd*(2*l5 + 2*s) + 2*l2*thetad*sin(phi - theta)*(l5 + s)))/((l5 + s)^2 + l2^2 - 2*l2*cos(phi - theta)*(l5 + s))^2;



end
