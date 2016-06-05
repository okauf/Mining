% Display motion for given control and configuration


%%%%%
% Set up fixed values of the system
%%%%%

initiateParameters;

%%
% ODE is to be solved on the intervall [0,T]
T = 1;
t_span = [0,T];

%%
% initial value of [s,theta,sd,thetad]
x_0 = [3; -pi/8; 0; 0];

%%
% torques defined as piecewise constant functions
% Later, in Discretization, use supporting points

u1 = [  -10000;
        -10000;
        -10000];
u2 = [  -10000;
        -10000;
        -10000];
tau_B1 = piecwConst(u1,T);
tau_B2 = piecwConst(u2,T);



%%%%%
% Set up discretization of the system
%%%%%

%%
% time points for discretization:
% t_pts = [0=t_0,t_1,...,t_N=T] = [0,T/N,2T/N,...,(N-1)T/N,T]
% mesh = [t_1-t_0,t_2-t_1,...,t_N-t_{N-1}] = [T/N,...,T/N]
N = 40;                  % number of time steps
t_pts = (0:N)*T/N;       % time points of discretization
mesh = ones(1,N)*T/N;    % time steps of discretization

%%
% Solve the ODE
% evaluate the solution on the given discretization
% x_opt is then the discrete state for the given parameters
% x_opt will be used for the optimization as reference state
sol_ode = solveLagrODE(params,x_0,t_span,tau_B1,tau_B2);
x_opt = deval(sol_ode,t_pts);

[lr1,lr2] = stheta2lr(params,x_opt(1,:),x_opt(2,:));

DataVisualization(params,lr1,lr2);
