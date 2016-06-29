function [ objFct, x_opt, trajForP ] = exampleConf4(params, N, rkOrder, ref)
% Example configuration 1
% Input:
%   params      optimal parameters
%   N           number of time steps in discretization
%   rkOrder     index of used Runge Kutta method
%               1: explicit Euler
%               2: RK of 2nd order
%               3: RK of 3nd order / Simpson's rule
%               4: RK of 4nd order
%   ref        index for reference for comparing x(p) with
%               1: x_opt
%               2: x(p_opt)
% Output:
%   objFct      objective function for this problem
%   x_opt       reference trajectory
%   trajForP    function to calculate a trajectory for given p
%
% The parameters come from initiateParameters and are thus fixed


%%%%%
% Set up fixed values of the system
%%%%%

% actual parameters, corresponding to the actual state
p_opt = pFromParams(params);
fix_params = nonOptParams(params);

%%
% ODE is to be solved on the intervall [0,T]
T = 10;
t_span = [0,T];

%%
% initial value of [s,theta,sd,thetad]
x_0 = [3; -pi/4; 0; 0];

%%
% torques defined as piecewise constant functions
% Later, in Discretization, use supporting points
u1 = [  -30000;
        -27000];
u2 = [  -9000;
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
t_pts = (0:N)*T/N;       % time points of discretization
mesh = ones(1,N)*T/N;    % time steps of discretization

%%
% Solve the ODE
% evaluate the solution on the given discretization
% x_opt is then the discrete state for the given parameters
% x_opt will be used for the optimization as reference state
sol_ode = solveLagrODE(params,x_0,t_span,tau_B1,tau_B2);
x_opt = deval(sol_ode,t_pts);

%%
% discrete torque values to avoid repetitive evaluation
tau_B1_val = tau_B1(t_pts);
tau_B2_val = tau_B2(t_pts);
u = [tau_B1_val(:)';
     tau_B2_val(:)'];


%%
% Approximation function
% calculate state approximation depending on the parameters
% choose from explicit Euler, RK order 2, 3 and 4 with index 1/2/3/4
approxFct = chooseApprox(fix_params,x_opt,u,mesh,rkOrder);

%%
% Objective function of tracking type
% Compare approximation from p to optimal approximation
% 
% |x_opt - x(p_opt)| gives an approximation error.
% To avoid this error in the objective function, one may
% use |x(p_opt) - x(p)| as the difference of the approximations
if ref == 1
    x_ref = x_opt;
else
    x_ref = approxFct(p_opt);
end
objFct = trackingTypeFct(x_ref, approxFct);


%%
% Trajectory for given p.
% As this type of problem only approximates p and does not
% explicitely create a trajectory purely dependent on p,
% one should also compare a resulting trajectory with the
% real one.
    function xp = trajForPFct(p)
        % Compute a trajectory for given p.
        pparams = pToParams(p,fix_params);
        sol_ode_p = solveLagrODE(pparams,x_0,t_span,tau_B1,tau_B2);
        xp = deval(sol_ode_p,t_pts);
    end

trajForP = @(p) trajForPFct(p);


end
