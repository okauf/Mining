function [ res_sqp, relerr, elapsed_time ] =  testNScale(params, filename)
% Compare solutions using different discretization N
% N ∊ {50,100,150,200,250,300,350,400,450,500}
% Test different N for different RK-methods
% The objective function also gets scaled with different K
% K ∊ {1,1000,10000,100000}
% One common starting point is chosen for all problems

p_opt = pFromParams(params);

% deviation
dev = 0.3;
p_0 = pToRandP(p_opt,dev);

exampleConf = @(a,b,c) exampleConf4(a,b,c,1); % always x_opt as reference
lb = zeros(size(p_0));
rkOrder = [1,2,3,4];
N = [50,200,500];
K = [1,1000];
descr = 'results res{rkOrder}{N}{K} contain p,fval,exitflag,output,rkOrder,N,K';

%%
% Setup for sqp solver
options                             = optimoptions('fmincon');
%options.Display                     = 'iter';
options.SpecifyObjectiveGradient    = true;
options.Algorithm                   = 'sqp';
options.OptimalityTolerance         = 1e-10;
options.StepTolerance               = 1e-10;
options.OutputFcn                   = @outfun;

% Every Runge Kutta method and every N gives a different output of the solver.
% Save it in an object res_sqp for later plots.

for n=1:length(N)
    for r=1:length(rkOrder)
        for k=1:length(K)
            disp(['status: n=', num2str(n), ' r=', num2str(r), ' k=', num2str(k)]);
            % create test instance
            [objFct, x_ref, trajForP] = exampleConf(params,N(n),rkOrder(r));
            objFct = rescaleObjFct(objFct,K(k));

            % solve and time
            tic;
            [p,fval,exitflag,output] = fmincon(objFct,p_0,[],[],[],[],lb,[],[], options);
            elapsed_time(n,r,k) = toc;

            % collect results
            res_sqp{r}{n}{k} = {p,fval,exitflag,output,rkOrder(r),N(n),K(k)};
        end
    end
end

if nargin >= 2
    try
        save(filename,'res_sqp','p_opt','p_0','rkOrder','N','K','dev','descr','elapsed_time','relerr');
    catch
        disp('failed to save results');
    end
end


% Function evaluated at every iteration
% collecting relative errors of the trajectory
    function stop = outfun(p,optimValues,state)
        stop = false;
        %relerr = norm(x_ref - trajForP(p))/norm(x_ref);
        %plot(optimValues.iteration,relerr);
        relerr(n,r,k,optimValues.iteration+1) = norm(x_ref - trajForP(p))/norm(x_ref);
    end



end
