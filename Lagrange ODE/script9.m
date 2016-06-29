clear;

initiateParameters;
p_opt = pFromParams(params);

N = 100;
rkOrder = 1;
ref = 1;
deviation = 0.5;
instance = [N;rkOrder;ref;deviation];

objFct = exampleConf4(params,N,rkOrder,ref);

lb = zeros(10,1);

options                 = optimoptions('fmincon');
options.Display         = 'iter';
options.SpecifyObjectiveGradient = true;
options.Algorithm       = 'sqp';
%options.Algorithm       = 'interior-point';
options.OptimalityTolerance = 1e-10;
options.StepTolerance   = 1e-10;
%options.CheckGradients  = true;
%options.ScaleProblem    = 'obj-and-constr';    % often not converging

for i=1:1
p_0 = pToRandP(p_opt,deviation);
[p,fval,exitflag,output] = fmincon(objFct,p_0,[],[],[],[],lb,[],[], options);
res(i,:) = [output.iterations;
            norm(p-p_opt)/norm(p_opt);
            exitflag]';
end

instance'
res
