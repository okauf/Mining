%function [p,p_opt,p_0,objFct,newObj] = script9(params)
p_opt = pFromParams(params);
%p_0 = ones(10,1);
p_0 = 0.9*p_opt;

objFct = exampleConf2(params,40,4);

lb = zeros(size(p_0));

options                 = optimoptions('fmincon');
options.Display         = 'iter';
options.SpecifyObjectiveGradient = true;
options.Algorithm       = 'sqp';
options.StepTolerance   = 1e-10;
options.CheckGradients  = true;
options.ScaleProblem    = 'obj-and-constr';


[p,fval,exitflag,output] = fmincon(objFct,p_0,[],[],[],[],lb,[],[], options);

%{
function [f,g] = scaledFct(p)
    [f,g] = objFct(100*p);
    g = g*100;
end
newObj = @(p) scaledFct(p);

p = fmincon(newObj,p_0/100,[],[],[],[],lb,[],[], options);
%}
%end
