initiateParameters;

N = 20;
[objFct,constrFct,x_opt] = multShootProblem(params,N);

options                 = optimoptions('fmincon');
options.Display         = 'iter';
options.SpecifyObjectiveGradient = true;
options.SpecifyConstraintGradient= true;
options.Algorithm       = 'sqp';
options.OptimalityTolerance = 1e-10;
options.StepTolerance   = 1e-10;
%options.CheckGradients  = true;

%p_0 = pToRandP(p_opt,deviation);
z_0 = [x_opt(:);10*ones(10,1)];
[z,fval,exitflag,output] = fmincon(objFct,z_0,[],[],[],[],[],[],constrFct,options);
