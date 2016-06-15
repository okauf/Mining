% Calculate symbolic derivatives of A and b from lagrdd

initiateParameters;

syms M1 M2 I_B1 I_B2 I_P1 I_P2 mu_B1 mu_B2 mu_P1 mu_P2;
p = [M1 M2 I_B1 I_B2 I_P1 I_P2 mu_B1 mu_B2 mu_P1 mu_P2];
params = pToParams(p,params);
syms l2 l3 r_B1 r_B2 r_P1 r_P2 l4 l5 r_M2;
params.r_B1 = r_B1;
params.r_B2 = r_B2;
params.r_P1 = r_P1;
params.r_P2 = r_P2;
params.l4 = l4;
params.l5 = l5;
params.r_M2 = r_M2;
params.l2 = l2;
params.l3 = l3;

syms a_s b_s c_s a_theta b_theta c_theta;
syms x1 x2 x3 x4;
x = [x1 x2 x3 x4];
%syms a_s(x1,x2,x3,x4) b_s(x1,x2,x3,x4) c_s(x1,x2,x3,x4) a_theta(x1,x2,x3,x4) b_theta(x1,x2,x3,x4) c_theta(x1,x2,x3,x4);

syms Q_s dTds dVds Q_theta dTdtheta dVdtheta;

%%
%%%%%%%%%%%%%%%%
% lagrdd
matA = [params.M1 + params.M2 + params.I_B2/params.r_B2^2 + params.I_P2/params.r_P2^2 + 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*a_s, ...
            0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*b_s;
            0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*a_theta, ...
            params.M1*(x(1) + params.l5)^2 + params.M2*(0.25*params.r_M2^2 + 1/12*params.l4^2 + 2*(x(1) + params.l5 - params.l4/2)^2) + 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*b_theta];


%%
% Q_s, dTds, ... also depend on x, but outsorced
bVec = [Q_s + dTds - dVds ...
            - 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*c_s;
        Q_theta + dTdtheta - dVdtheta ...
            - (2*params.M1*(x(1) + params.l5) + 4*params.M2*( x(1) + params.l5 - 0.5*params.l4))*x(3)*x(4) ...
            - 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*c_theta];

%%
%%%%%%%%%%%%%%%%
% Energies
syms dlr1d2ds(s,theta,sd,thetad) dlr1d2dtheta(s,theta,sd,thetad) s theta sd thetad;
syms g;
params.g = g;
x2 = [s theta sd thetad];

%Derivative of T wrt s
dTds = (params.M1*(s+params.l5) + 2*params.M2*(s+params.l5 - ...
    0.5*params.l4))* thetad^2 + 0.5*( params.I_B1/params.r_B1^2 + ...
    params.I_P1/params.r_P1^2) * dlr1d2ds;

%Derivative of V wrt s
dVds = (params.M1+params.M2)*params.g*sin(theta);

%Derivative of T wrt theta
dTdtheta = 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2) * dlr1d2dtheta;

%Derivative of V wrt theta
dVdtheta = (params.M1 *(s+params.l5) + params.M2 * (s+params.l5 ...
    - params.l4/2))*params.g*cos(theta);



%%
%%%%%%%%%%%%%%%%
% Forces
syms phi tau_B1 tau_B2
lr1d = ((s + params.l5)*sd - params.l2*(sd*cos(phi-theta) + (s + params.l5)*sin(phi-theta)*thetad)) / ...
    sqrt(params.l2^2 + (s + params.l5)^2 - 2*params.l2*(s + params.l5)*cos(phi-theta));

r1 = (s+params.l5)*[cos(theta);sin(theta)] - params.l2*[cos(phi);sin(phi)];
r1 = r1 / norm(r1);
r2 = [cos(theta);sin(theta)];


% rA1 is the position vector of attachement point A1 for force FA1
% drA1ds and drA1dtheta are the derivatives w.r.t. s and theta
drA1ds = [cos(theta);sin(theta)];
drA1dtheta = (s + params.l5)*[-sin(theta);cos(theta)];

% rA2 is the position vector of attachement point A2 for force FA2
% drA2ds and drA2dtheta are the derivatives w.r.t. s and theta
drA2ds = [cos(theta);sin(theta)];
drA2dtheta = s*[-sin(theta);cos(theta)];


% forces FA1 and FA2
FA1 = (tau_B1/params.r_B1 - params.mu_B1 * lr1d/params.r_B1 - params.mu_P1 * lr1d/params.r_P1) * r1;
FA2 = (tau_B2/params.r_B2 - params.mu_B2 * sd/params.r_B2 - params.mu_P2 * sd/params.r_P2) * r2;

% generalized forces
Q_s = drA1ds'*FA1 + drA2ds'*FA2;
Q_theta = drA1dtheta'*FA1 + drA2dtheta'*FA2;


%%
%%%%%%%%%%%%%%%%
% derRopeLength

% Z is the numerator of lr1d2
Z = (s+params.l5)^2 * sd^2 - 2*(s+params.l5) * sd^2 * params.l2*cos(phi-theta) ...
    - 2*(s+params.l5)^2 *sd * thetad * params.l2 * sin(phi - theta) ...
    + params.l2^2 * sd^2 * cos(phi-theta)^2 + 2*params.l2^2*(s+params.l5) * sd * thetad * cos(phi-theta) * sin(phi-theta) ...
    + params.l2^2 * (s+params.l5)^2 * thetad^2 * sin(phi-theta)^2;

% N is the denominator of lr1de
N = params.l2^2 + (s+params.l5)^2 - 2*params.l2*(s+params.l5) * cos(phi-theta);

%Derivative of dot(lr1)^2 wrt s
dlr1d2ds = ((2*(s+params.l5)*sd^2 - 2*sd^2*params.l2 * cos(phi-theta) - ...
    4*(s+params.l5) * sd * thetad * params.l2 * sin(phi-theta) ...
    + 2*params.l2^2 * sd * thetad * cos(phi-theta) * sin(phi-theta) + 2 * ...
    params.l2^2 * (s+params.l5) * thetad^2 * sin(phi-theta)^2) * N ...
    - Z * (2*(s+params.l5) - 2*params.l2 * cos(phi-theta))) / N^2;

%Derivative of dot(lr1)^2 wrt theta
dlr1d2dtheta = ((-2*(s+params.l5) * sd^2 * params.l2 * sin(phi-theta) + 2*(s+params.l5)^2 * ...
    sd * thetad * params.l2 * cos(phi-theta) + 2*params.l2^2 * sd^2 * cos(phi-theta) * sin(phi-theta) ... 
    + 2*params.l2^2*(s+params.l5) * sd * thetad * sin(phi-theta)^2 ...
    - 2*params.l2^2*(s+params.l5) * sd * thetad * cos(phi-theta)^2 ...
    - 2*params.l2^2*(s+params.l5)^2 * thetad^2 * sin(phi-theta) * cos(phi-theta)) * N ...
    - Z * (-2*params.l2*(s+params.l5) * sin(phi-theta))) / N^2;



%%
%%%%%%%%%%%%%%%
% derderRopeLength

N = params.l2^2 + (s+params.l5)^2 - 2*params.l2*(s+params.l5) * cos(phi-theta);

a_s = (2*(s+params.l5)^2 - 4*(s+params.l5)*params.l2 * cos(phi-theta) + 2*params.l2^2 * cos(phi-theta)^2) / N;

b_s = (-2*(s+params.l5)^2*params.l2 * sin(phi-theta) + 2*params.l2^2*(s+params.l5) * cos(phi-theta) * sin(phi-theta))/N;

c_s = ((4* sd^2*(s+params.l5) - 4*sd^2*params.l2 * cos(phi-theta) - ...
    4*sd*thetad*(s+params.l5)*params.l2 * sin(phi-theta) ...
    - 4*sd*thetad*(s+params.l5)*params.l2 * sin(phi-theta) ...
    + 2*thetad^2*(s+params.l5)^2 * params.l2 * cos(phi-theta) ...
    + 6*sd*thetad*params.l2^2 * cos(phi-theta) * sin(phi-theta) ...
    + 2*thetad^2*params.l2^2*(s+params.l5) * sin(phi-theta)^2 - 2*thetad^2 * params.l2^2 * ... 
    (s+params.l5) * cos(phi-theta)^2) * N ...
    - (2*sd*(s+params.l5)^2 - 4*sd*(s+params.l5)*params.l2 * cos(phi-theta) ...
    - 2 * thetad*(s+params.l5)^2*params.l2 * sin(phi-theta) ...
    + 2*sd*params.l2^2 * cos(phi-theta)^2 + 2*thetad*params.l2^2*(s+params.l5) * ... 
    cos(phi-theta) * sin(phi-theta))* (2*sd*(s+params.l5) - 2*sd*params.l2 * cos(phi-theta) -2*thetad * ...
    params.l2*(s+params.l5) * sin(phi-theta))) / N^2;



a_theta = (-2*(s+params.l5)^2 * sin(phi-theta) + 2*params.l2^2*(s+params.l5) * ...
    cos(phi-theta) * sin(phi-theta)) / N;

b_theta = 2*params.l2^2*(s+params.l5)^2 * sin(phi-theta)^2 / N;

c_theta = ((-4*(s+params.l5)*sd^2*params.l2 * sin(phi-theta) ...
    + s*(s+params.l5)^2*sd*thetad*params.l2 * cos(phi-theta) ...
    + 2*params.l2^2 * sd^2 * cos(phi-theta) * sin(phi-theta) ...
    + 6*params.l2^2 * (s+params.l5)*sd*thetad * sin(phi-theta)^2 ...
    - 2*params.l2^2 * (s+params.l5)*sd*thetad * cos(phi-theta)^2 ...
    - 2*params.l2^2 * (s+params.l5)^2 * thetad^2 * sin(phi-theta) * cos(phi-theta)) * N ...
    - (-2*(s+params.l5)^2 *sd *params.l2 * sin(phi-theta) + 2*params.l2^2*(s+params.l5)*sd * ...
    cos(phi-theta) * sin(phi-theta) + 2*params.l2^2*(s+params.l5)^2 * thetad * sin(phi-theta)^2) * ...
    (2*(s+params.l5)* sd - 2*params.l2*sd * cos(phi-theta) ...
    - 2*params.l2*thetad*(s+params.l5) * sin(phi-theta))) / N^2;





for i = 1:4
    x2(i)
    simplify(diff(c_theta,x2(i)))
    pause;
end
