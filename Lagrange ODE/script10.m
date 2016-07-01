initiateParameters;
load('data/compareErr35_1.mat');

N = 1000;
[objFct,constrFct,x_opt,z_0] = multShootProblem(params,p_0,N);
%objFct = rescaleObjFct(objFct,1/1000);

disp(['came so far']);

options                 = optimoptions('fmincon');
options.Display         = 'iter';
options.SpecifyObjectiveGradient = true;
options.SpecifyConstraintGradient= true;
%options.Algorithm       = 'sqp';
options.Algorithm       = 'interior-point';
%options.Algorithm       = 'active-set';
options.OptimalityTolerance = 1e-10;
options.StepTolerance   = 1e-10;
options.FunctionTolerance   = 1e-10;
%options.CheckGradients  = true;
options.OutputFcn       = @optimplotfval;
options.UseParallel     = true;

%p_0 = pToRandP(p_opt,deviation);
%z_0 = [x_opt(:);10*ones(10,1)];
disp(['starting']);
[z,fval,exitflag,output] = fmincon(objFct,z_0,[],[],[],[],[],[],constrFct,options);
