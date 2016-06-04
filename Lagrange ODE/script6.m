% Create instance for optimization


%%%%%
% Set up the Optimization Problem
%%%%%

initiateParameters;
p_opt = pFromParams(params);


%%
% Choose initial values for the parameters
% They differ from the real, unknown parameters
% p_0 = [M1,M2,I_B1,I_B2,I_P1,I_P2,mu_B1,mu_B2,mu_P1,mu_P2]
p_0 = 0.5*p_opt;


%%
% Create an objective function from several instances
% Every instance has a state to approximate,
% obj = 0.5*|x(p_opt) - x(p)|^2
% Here, optimize p such that all of them are a close fit
obj1 = exampleConf1(params,100,1);
obj2 = exampleConf2(params,100,4);
obj3 = exampleConf3(params,100,3);

objFct = @(p) obj1(p) + obj2(p) + obj3(p);


%%
% Constraint
% Non-negativity lower bound for parameters
lb = zeros(size(p_0));



%%%%%
% Solve the Optimization Problem
%%%%%


%%
tic;
options                 = optimoptions('fmincon');
options.Display         = 'iter';
% options.GradObj         = 'on';
% options.GradConstr      = 'on';
% options.Hessian         = 'user-supplied';
% options.HessFcn         = @hessianMap;

p = fmincon(objFct,p_0,[],[],[],[],lb,[],[], options);
toc;
