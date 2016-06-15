% Create instance for optimization

%%%%%
% Set up fixed values of the system
%%%%%

initiateParameters;

% actual parameters, corresponding to the actual state
p_opt = pFromParams(params);

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
u1 = [  -20000;
        -15000;
        -50000;
        -30000;
        -35000;
        -40000;
        -10000;
        -50000;
        -40000;
        -30000;
        -10000;
        -10000];
u2 = [  -15000;
        -10000;
        -8000;
        0;
        0;
        -10000;
        -20000;
        -10000;
        -5000;
        -5000;
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
N = 100;                  % number of time steps
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
% choose from explicit Euler, RK order 2, 3 and 4 with index 1/2/3/4
approxFct = chooseApprox(fix_params,x_opt,u,mesh,1);

%%
% Objective function of tracking type
% Compare approximation from p to optimal approximation
% 
% |x_opt - x(p_opt)| gives an approximation error.
% To avoid this error in the objective function, 
% use |x(p_opt) - x(p)| as the difference of the approximations
objFct = trackingTypeFct(approxFct(p_opt), approxFct);

%%
% Constraint
% Non-negativity lower bound for parameters
lb = zeros(size(p_0));

%{
%%
% Gradient of objective functin w.r.t. p
% Calculate derivatives w.r.t. p
% Do this by symbolically deriving the objective function

syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10;
p_sym = [ p1, p2, p3, p4, p5, p6, p7, p8, p9, p10];

% approximation function useable for symbolics
approxFctSym = explEulerApproxSym(fix_params,x_opt,u,mesh);

% objective function useable for symbolics
objFctSym = trackingTypeFct(approxFctSym(p_opt), approxFctSym);

% evaluate objective function at symbolic p
y_sym = objFctSym(p_sym);
y_sym = simplify(y_sym);

% calculate derivative w.r.t. every symbolic parameter
for i=1:length(p_sym)
    yd_sym(i) = diff(y_sym,p_sym(i));
end

yd_sym = simplify(yd_sym);

% Gradient function of the objective
objGrad = @(p) double(subs(yd_sym(:),p_sym(:),p(:)));

% For optimizatoin, objective and gradient
% It evaluates p and returns the objective function,
% and if requested also return the gradient.
objWithGrad = getObjWithGrad(objFct,objGrad);


%%%%%
% Solve the Optimization Problem
%%%%%


%%
tic;
options                 = optimoptions('fmincon');
options.Display         = 'iter';
options.SpecifyObjectiveGradient = true;
% options.GradObj         = 'on';
% options.GradConstr      = 'on';
% options.Hessian         = 'user-supplied';
% options.HessFcn         = @hessianMap;

p = fmincon(objWithGrad,p_0,[],[],[],[],lb,[],[], options);
toc;
%}
