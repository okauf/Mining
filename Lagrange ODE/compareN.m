function compareN(params)
% Compare solutions using different discretization N
% N ∊ {5,10,20,50,100}
% Test different N for different RK-methods, using different solvers

p_opt = pFromParams(params);
p_0 = 0.5*p_opt;
exampleConf = @(a,b,c) exampleConf2(a,b,c);
lb = zeros(size(p_0));
rkOrder = [1,2,3,4];
N = [10,20,50,100];
descr = 'results res{rkOrder}{N}';

%%
% Setup for sqp solver
options                 = optimoptions('fmincon');
options.Display         = 'iter';
options.SpecifyObjectiveGradient = true;
options.Algorithm       = 'sqp';

% Every Runge Kutta method and every N gives a different output of the solver.
% Save it in an object res_sqp for later plots.

for k=1:length(rkOrder)
    for n=1:length(N)
        objFct = exampleConf(params,N(n),rkOrder(k));
        [p,fval,exitflag,output] = fmincon(objFct,p_0,[],[],[],[],lb,[],[], options);
        res_sqp{k}{n} = {p,fval,exitflag,output,rkOrder(k),N(n)};
    end
end


%%
% Setup for interior-point solver
options                 = optimoptions('fmincon');
options.Display         = 'iter';
options.SpecifyObjectiveGradient = true;
options.Algorithm       = 'interior-point';

% Every Runge Kutta method and every N gives a different output of the solver.
% Save it in an object res_int for later plots.

for k=1:length(rkOrder)
    for n=1:length(N)
        objFct = exampleConf(params,N(n),rkOrder(k));
        [p,fval,exitflag,output] = fmincon(objFct,p_0,[],[],[],[],lb,[],[], options);
        res_int{k}{n} = {p,fval,exitflag,output,rkOrder(k),N(n)};
    end
end

save('data/compareN.mat','res_sqp','res_int','p_opt','p_0','rkOrder','N','descr');


end
