initiateParameters;

%tau_p1  = @(t) params.I_cr*sin(t)-0.1*t;
%tau_p2  = @(t) -params.I_cr*sin(t);
%tau_p1  = @(t) params.I_cr*0.1*t;
%tau_p2 = @(t) -params.I_cr*0.1*sqrt(t);

tau_p1 = @(t) params.I_cr*( [t<2*pi]*sin(t) - [t>=5 && t<8]*(t-5) + [t>=8 && t<10]*1 + [t>=10]*2 );
tau_p2 = @(t) params.I_cr*( -[t<pi]*sin(t) + [t>=pi && t<2*pi]*0.1*sin(t) + [t>=2*pi] );

[lr1,lr2,sol] = strTwo(params,tau_p1,tau_p2,[0,14]);
A = DataVisualization(params,lr1,lr2);
