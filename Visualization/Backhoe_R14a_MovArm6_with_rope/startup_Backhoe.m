% Copyright 2012-2014 The MathWorks, Inc.
Backhoe_Arm_HomeDir = pwd;

addpath(pwd);
addpath([pwd '\Scripts']);
addpath([pwd '\Libraries']);
addpath([pwd '\Images']);
addpath([pwd '\Images\Diagrams']);
addpath([pwd '\Actuation']);
addpath([pwd '\Custom_Orifice']);
addpath([pwd '\Custom_Valve']);
addpath([pwd '\Param_Est']);
addpath([pwd '\Optimize']);

if(exist('Custom_Orifice')==7)
    cd Custom_Orifice
    if((exist('+Hydraulic')==7) && ~exist('Hydraulic_lib'))
        disp('Building Custom Simscape Library...');
        ssc_build Hydraulic
        disp('Finished Building Library.');
    end
    cd ..
end

Backhoe_Model_PARAM
Backhoe_Arm
open('Backhoe_Demo_Script.html')


t = 1:20; t = t';
s = zeros(size(t));
sd = diff(s); sd = [sd; sd(end)];
sdd = diff(sd); sdd = [sdd; sdd(end)];

theta = linspace(0,pi/4,20)';
thetad = diff(theta); thetad = [thetad; thetad(end)];
thetadd = diff(thetad); thetadd = [thetadd; thetadd(end)];

a = 5 ;%+ 0.45; %0.39;
c = 5.625 - s; %(0.5 + 4.1 + 0.5 + 0.2) 4.8
beta = pi/2 + theta;

b = sqrt(-2*a*c.*cos(beta) + a^2 + c.^2);
ang = acos((-a*cos(beta) + c)./b) + 3/2*pi;

angd = diff(ang); angd = [angd; angd(end)];
angdd = diff(angd); angdd = [angdd; angdd(end)];

