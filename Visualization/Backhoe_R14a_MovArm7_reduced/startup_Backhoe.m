% Copyright 2012-2014 The MathWorks, Inc.
Backhoe_Arm_HomeDir = pwd;

addpath(pwd);
addpath([pwd '\Scripts']);
addpath([pwd '\Libraries']);

Backhoe_Model_PARAM
Backhoe_Arm

t = 1:20; t = t';
s = zeros(size(t));
sd = diff(s); sd = [sd; sd(end)];
sdd = diff(sd); sdd = [sdd; sdd(end)];

theta = linspace(0,pi/4,20)';
thetad = diff(theta); thetad = [thetad; thetad(end)];
thetadd = diff(thetad); thetadd = [thetadd; thetadd(end)];

a = 5 ;
c = 5.625 - s;
beta = pi/2 + theta;

b = sqrt(-2*a*c.*cos(beta) + a^2 + c.^2);
ang = acos((-a*cos(beta) + c)./b) + 3/2*pi;

angd = diff(ang); angd = [angd; angd(end)];
angdd = diff(angd); angdd = [angdd; angdd(end)];