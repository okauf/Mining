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
x_0 = [3; -pi/4; 0; 0];

%%
% torques defined as piecewise constant functions
% Later, in Discretization, use supporting points
u1 = [-25000];
u2 = [-25000];
tau_B1 = piecwConst(u1,T);
tau_B2 = piecwConst(u2,T);



%%%%%
% Set up discretization of the system
%%%%%

%%
% time points for discretization:
% t_pts = [0=t_0,t_1,...,t_N=T] = [0,T/N,2T/N,...,(N-1)T/N,T]
% mesh = [t_1-t_0,t_2-t_1,...,t_N-t_{N-1}] = [T/N,...,T/N]
N = 20;                  % number of time steps
t_pts = (0:N)*T/N;       % time points of discretization
mesh = ones(1,N)*T/N;    % time steps of discretization

%%
% Solve the ODE
% evaluate the solution on the given discretization
% x_opt is then the discrete state for the given parameters
% x_opt will be used for the optimization as reference state
sol_ode = solveLagrODE(params,x_0,t_span,tau_B1,tau_B2);
x_opt = deval(sol_ode,t_pts);   %%%what ordering?

%%
% discrete torque values to avoid repetitive evaluation
tau_B1_val = tau_B1(t_pts);
tau_B2_val = tau_B2(t_pts);
u = [tau_B1_val(:)';
     tau_B2_val(:)'];



%%%%%
% Set up the Optimization Problem
%%%%%

%%
% Choose initial values for the parameters
% They differ from the real, unknown parameters
% p_0 = [M1,M2,I_B1,I_B2,I_P1,I_P2,mu_B1,mu_B2,mu_P1,mu_P2]
p_0 = ones(10,1);   %%%
fix_params = nonOptParams(params);

%%
% Approximation function
% calculate state approximation depending on the parameters
approxFct = explEulerApprox(fix_params,x_opt,u,mesh);

%%
% Objective function of tracking type
% difference between x_opt and discrete approx from parameters
objFct = trackingTypeFct(x_opt, approxFct);



%%%%%
% Solve the Optimization Problem
%%%%%

%{
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
%}
