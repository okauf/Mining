initiateParameters;

tau_p1  = @(t) params.I_cr*sin(t)-0.1*t;
tau_p2  = @(t) -params.I_cr*sin(t);
%tau_p1  = @(t) params.I_cr*0.1*t;
%tau_p2 = @(t) -params.I_cr*0.1*sqrt(t);

[lr1,lr2,sol] = strTwo(params,tau_p1,tau_p2);
DataVisualization(params,lr1,lr2);
