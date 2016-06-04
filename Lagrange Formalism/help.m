function [ Der ] = help( )

syms s sd theta thetad sdd thetadd l5 l2 phi

Z = (2*(s-l5)*sd - 2*l2*(sd*cos(phi-theta)+(s+l5)*sin(phi-theta)*thetad))*l2*(s+l5)*sin(phi-theta);

N = l2^2 + (s+l5)^2 - 2*l2*(s+l5)*cos(phi-theta);

dZ = (2*sd^2+2*(s+l5)*sdd-2*l2*(sdd*cos(phi-theta)+sd*sin(phi-theta)*thetad + sd*sin(phi-theta)*thetad - (s+l5)*cos(phi-theta)*thetad^2 + (s+l5)*sin(phi-theta)*thetadd)) * ...
    l2*(s+l5)*sin(phi-theta)...
    +(2*(s+l5)*sd - 2*l2*(sd*cos(phi-theta)+(s+l5)*sin(phi-theta)*thetad)) * ...
    (l2*sd*sin(phi-theta)-l2*(s+l5)*cos(phi-theta)*thetad);

dN = 2*(s+l5)*sd - 2*l2*sd*cos(phi-theta)-2*l2*(s-l5)*sin(phi-theta)*thetad;

Der = (N*dZ-Z*dN)/N^2;

Der = simplify(Der);

end

