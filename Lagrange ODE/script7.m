% Calculate symbolic derivatives of A and b from lagrdd

initiateParameters;

syms M1 M2 I_B1 I_B2 I_P1 I_P2 mu_B1 mu_B2 mu_P1 mu_P2;
p = [M1 M2 I_B1 I_B2 I_P1 I_P2 mu_B1 mu_B2 mu_P1 mu_P2];
params = pToParams(p,params);
syms r_B1 r_B2 r_P1 r_P2 l4 l5 r_M2;
params.r_B1 = r_B1;
params.r_B2 = r_B2;
params.r_P1 = r_P1;
params.r_P2 = r_P2;
params.l4 = l4;
params.l5 = l5;
params.r_M2 = r_M2;


syms a_s b_s c_s a_theta b_theta c_theta;
syms x1 x2 x3 x4;
x = [x1 x2 x3 x4];

syms Q_s dTds dVds Q_theta dTdtheta dVdtheta;

matA = [params.M1 + params.M2 + params.I_B2/params.r_B2^2 + params.I_P2/params.r_P2^2 + 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*a_s, ...
            0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*b_s;
            0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*a_theta, ...
            params.M1*(x(1) + params.l5)^2 + params.M2*(0.25*params.r_M2^2 + 1/12*params.l4^2 + 2*(x(1) + params.l5 - params.l4/2)^2) + 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*b_theta];


%{
for k = 1:10
    p(k)
    diff(matA,p(k))
    pause;
end
%}

bVec = [Q_s + dTds - dVds ...
            - 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*c_s;
        Q_theta + dTdtheta - dVdtheta ...
            - (2*params.M1*(x(1) + params.l5) + 4*params.M2*( x(1) + params.l5 - 0.5*params.l4))*x(3)*x(4) ...
            - 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*c_theta];

%{
% only displays terms not related to T,V and Q
for k = 1:10
    p(k)
    diff(bVec,p(k))
    pause;
end
%}


syms dlr1d2ds dlr1d2dtheta s theta sd thetad;
syms g;
params.g = g;

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



for k = 1:10
    p(k)
    diff(dVdtheta,p(k))
    pause;
end

