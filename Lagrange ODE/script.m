initiateParameters;

%%
% grid for discretization
N = 20;                  % number of supporting points
int = 4;                 % length time interval
mesh = ones(1,N)*int/N;  % mesh used for discretization

%%
% initial value of [s,theta,sd,thetad]
x_0 = [3; -pi/4; 0; 0];

%%
% torques
k = 1500;   % fixed, such that the system barely moves
%tau_P1 = @(t) -k;
%tau_P2 = @(t) -k;
%tau_P1 = @(t) -k*([t<0.2*pi]*sin(t) - [t>=0.5 && t<0.8]*(t-0.5) + [t>=0.8]*1);
%tau_P2 = @(t) -k*([t<0.2*pi]*sin(t) - [t>=0.5 && t<0.8]*(t-0.5) + [t>=0.8]*1);
tau_P1 = @(t) -k - 30*k*t;
tau_P2 = @(t) -k - k*t;

%%
% state x
% x = [s; theta; sd; thetad]
% sd = ds/dt, thetad = d theta/dt
%
% control u
% u = [tau_P1; tau_P2]

% supporting points of state function
x0 = zeros(4,N+1); x = x0;

% supporting points of control function
u0 = zeros(2,N+1);

%%
tic;
options                 = optimoptions('fmincon');
options.Display         = 'iter';
% options.GradObj         = 'on';
% options.GradConstr      = 'on';
% options.Hessian         = 'user-supplied';
% options.HessFcn         = @hessianMap;

[ Cost, Constr ] = OptimizationFunctions(x, u, mesh);

x = fmincon(Cost,x0,[],[],[],[],[],[],Constr, options);
toc;

