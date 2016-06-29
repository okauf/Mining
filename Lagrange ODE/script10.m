initiateParameters;

N = 4;
[objFct,constrFct,x_opt] = multShootProblem(params,N);
objFct = rescaleObjFct(objFct,1/1000);

options                 = optimoptions('fmincon');
options.Display         = 'iter';
options.SpecifyObjectiveGradient = true;
%options.SpecifyConstraintGradient= true;
options.Algorithm       = 'sqp';
%options.Algorithm       = 'interior-point';
%options.Algorithm       = 'active-set';
options.OptimalityTolerance = 1e-10;
options.StepTolerance   = 1e-10;
options.FunctionTolerance   = 1e-10;
%options.CheckGradients  = true;

%p_0 = pToRandP(p_opt,deviation);
z_0 = [x_opt(:);10*ones(10,1)];
[z,fval,exitflag,output] = fmincon(objFct,z_0,[],[],[],[],[],[],constrFct,options);
