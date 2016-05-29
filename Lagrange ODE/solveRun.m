%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set initial values for ODE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initiateParameters;

% initial value of [s,theta,sd,thetad]
x_0 = [2; -pi/4; 1; 0.5];

% timespan for the ode
t_span = [0,1];

% torques
k = 25000;
tau_P1 = @(t) -k;
tau_P2 = @(t) -k;


%%%%%%%%%%%%%
% solve ODE %
%%%%%%%%%%%%%

sol = solveLagrODE(params, x_0, t_span, tau_P1, tau_P2);


%%%%%%%%%%%%%%%%%%%%%%%%
% translate to lr1,lr2 %
%%%%%%%%%%%%%%%%%%%%%%%%

s = sol.y(1,:);
theta = sol.y(2,:);

% calculate lr1-l3 via cosinus relation in general triangle
% c^2 := (lr1-l3)^2 = l2^2 + (s+l5)^2 - 2l2(s+l5)^2 * cos(phi-theta)
c_qu = params.l2^2 + (s+params.l5).^2 - 2*params.l2.*(s+params.l5).*cos(params.ang_base - theta);
c = sqrt(c_qu);
lr1 = c + params.l3;
lr2 = s + params.l1;


%%%%%%%%%%%%%
% visualize %
%%%%%%%%%%%%%

% function is in directory above
addpath ../

DataVisualization(params,lr1,lr2);
