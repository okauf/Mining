function [sol] = unstrTwoModTorques()

% for the computation of the angular acc of the cable reel at the base
% level the torque is composite of the control torque and the torque acting
% on the cable reel due to the force on the pulley (gravity)

% here, the length of the second rope is assumed to be constant
% therefore the decomposition of the force is only dependent on the length
% of the first rope

% solution for rope length and angular velocity takes imaginary values ->
% reasonable input function for torque control?

initiateParameters;

r_cr = 0.20;   % radius of cable reel
h_cr = 0.15;    % height of cable reel
rho_cr = 7860;  % [kg/m³] density of steel for cable reel
V_cr = pi*r_cr^2*h_cr;  % volume of cable reel
m_cr = rho_cr*V_cr;     % mass of cable reel
I_cr = 1/12*m_cr*(3*r_cr^2+h_cr^2);       % inertia of cable reel

lr2 = params.lr20;

alpha1 = @(lr1,lr2) acos(-((lr2-params.l1+params.l5)^2-params.l2^2-(lr1-params.l3)^2)/(2*params.l2*(lr1-params.l3)));
theta2 = @(lr1,lr2) acos(-((lr1-params.l3)^2-params.l2^2-(lr2-params.l1+params.l5)^2)/(2*params.l2*(lr2-params.l1+params.l5)));
beta_sum = @(lr1,lr2) pi-alpha1(lr1,lr2)-theta2(lr1,lr2);
%distinction for decomposition of forces!
beta1 = @(lr1,lr2) pi/2 - params.ang_base - alpha1(lr1,lr2);
beta2 = @(lr1,lr2) beta_sum(lr1,lr2) - beta1(lr1,lr2);
%tensile force on rope 1
F1 = @(lr1,lr2) params.F*sin(beta2(lr1,lr2))/sin(pi-beta1(lr1,lr2)-beta2(lr1,lr2));

% torque of pulley 1
tau_p1  = @(t) -2.9500e+03; %I_cr*10;%sin(t); % Input Function for torque on the cable reel at the base level

torque_r1 = @(x) r_cr * F1(x(1),lr2);%0.5 * 9.81 * 1000; %
rhs = @(t,x,u) [r_cr*x(2); (u+torque_r1(x(1)))/I_cr]; % ODE for cable 1
t_int = [0 1];      % time interval
x_0 = [23 0];       % start values of state variables

sol = ode45( @(t,x) rhs(t,x,tau_p1(t)), t_int, x_0);



end
