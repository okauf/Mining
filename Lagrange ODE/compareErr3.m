function [ absErr, relErr, absSTErr, relSTErr, fctVal, p_0, maxsdiff, maxtdiff ] =  compareErr3(params, filename)
% Compare convergence for different RK orders and different solvers
% functions to converge:
%   objective function
%   trajectory comparison, absolute and relative
%   trajectory comparison in s and θ, max absolute difference
% Choose common starting point for all problems.
% Fix discretization N=1000

p_opt = pFromParams(params);

% deviation
dev = 0.1;

exampleConf = @(a,b,c) exampleConf5(a,b,c,1); % always x_opt as reference
lb = zeros(10,1);
rkOrder = [1];
N = 1500;       % number of discretization time steps
K = 1000;       % scaling of the objective function
descr = 'result{alg,rkOrder,iter}';
algorithm = {'interior-point','sqp'};

%%
% Setup of the solvers, algorithm is set later
options                             = optimoptions('fmincon');
options.Display                     = 'iter';
options.SpecifyObjectiveGradient    = true;
options.OptimalityTolerance         = 1e-14;
options.StepTolerance               = 1e-14;
options.OutputFcn                   = @outfun;



p_0 = pToRandP(p_opt,dev);

for a=1:length(algorithm)
    options.Algorithm   = algorithm{a};
    for r=1:length(rkOrder)
        disp(['status: alg=', algorithm(a), ' r=', num2str(r)]);

        % create test instance
        [objFct, x_ref, trajForP] = exampleConf(params,N,rkOrder(r));
        objFct1 = exampleConf1(params,N,rkOrder(r),1);
        objFct2 = exampleConf2(params,N,rkOrder(r),1);
        objFct3 = exampleConf3(params,N,rkOrder(r),1);
        objFct4 = exampleConf4(params,N,rkOrder(r),1);
        combinedObjFct = combineObjFcts({objFct,objFct1,objFct2,objFct3,objFct4});
        scaledObjFct = rescaleObjFct(combinedObjFct,K);
        %scaledObjFct = rescaleObjFct(objFct,K);

        % solve
        tic;
        [p,fval,exitflag,output] = fmincon(scaledObjFct,p_0,[],[],[],[],lb,[],[], options);
        toc;

        % the results from the last iteration are not yet added to the results
        % do it manually
        optimValues.iteration = output.iterations;  % only value needed
        state = 0;                                  % not needed
        outfun(p,optimValues,state);
    end
end

if nargin >= 2
    try
        % also save the starting trajectory, optimized trajectory and reference trajectory
        % save in three separate files each corresponding s,sd,sdd, θ,θd,θdd
        % First: x_comp
        [x_comp,s,sd,sdd,theta,thetad,thetadd,t]  = trajForP(p);
        save('data/traj_x_comp.mat','x_comp','p','s','sd','sdd','theta','thetad','thetadd','t');

        % Second: x_p0
        [x_p0,s,sd,sdd,theta,thetad,thetadd,t]    = trajForP(p_0);
        save('data/traj_x_p0.mat','x_p0','p_0','s','sd','sdd','theta','thetad','thetadd','t');


        % Third: x_ref
        [x_ref,s,sd,sdd,theta,thetad,thetadd,t]    = trajForP(p_opt);
        save('data/traj_x_ref.mat','x_ref','p_opt','s','sd','sdd','theta','thetad','thetadd','t');


        save(filename,'absErr','relErr','absSTErr','relSTErr','fctVal','p_0','p','maxsdiff','maxtdiff');
    catch
        disp('failed to save results');
    end
end


% Function evaluated at every iteration
% Collect absolute and relative error of the full trajectory,
% and also of only the s and theta part of the trajectory.
% Also collect the value of the (unscaled) objective function.
    function stop = outfun(p,optimValues,state)
        if toc > 1800   % 30 min
            disp(['exceeded time limit']);
            stop = true;
        else
            stop = false;
        end

        iter = optimValues.iteration+1;
        trajP = trajForP(p);

        % x_ref and trajP are in R^4x(N+1)
        % Convert x_ref and trajP into vectors of the form
        % [s1,s2,...,θ1,θ2,...,sd1,sd2,...,θd1,θd2,...]
        x_ref_v = x_ref';
        x_ref_v = x_ref_v(:);
        trajP_v = trajP';
        trajP_v = trajP_v(:);

        % x_ref_st_v and trajP_st_v only contain s and θ part
        % [s1,s2,s3,...,θ1,θ2,θ3,...]
        x_ref_st_v = x_ref_v(1:2*(N+1));
        trajP_st_v = trajP_v(1:2*(N+1));

        absErr(a,r,iter)    = norm(x_ref_v - trajP_v,2);
        relErr(a,r,iter)    = norm(x_ref_v - trajP_v,2)/norm(x_ref_v,2);
        absSTErr(a,r,iter)  = norm(x_ref_st_v - trajP_st_v,2);
        relSTErr(a,r,iter)  = norm(x_ref_st_v - trajP_st_v,2)/norm(x_ref_st_v);
        fctVal(a,r,iter)    = objFct(p);
        maxsdiff(a,r,iter)  = max(abs(x_ref_v(1:(N+1)) - trajP_v(1:(N+1))));
        maxtdiff(a,r,iter)  = max(abs(x_ref_v((N+2):2*(N+1)) - trajP_v((N+2):2*(N+1))));
        disp(num2str(maxsdiff(a,r,iter)));
    end



end
