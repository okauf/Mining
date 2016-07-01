function [ objFct, constrFct, x_opt_vect, z_0 ] = multShootProblem(params, p_0, N)
% Multiple Shooting example problem setting
% Input:
%   params      optimal parameters
%   N           number of time steps in discretization
% Output:
%   objFct      objective function for this problem
%   constrFct   multiple shooting constraint
%   x_opt_vect  reference trajectory
%
% The parameters come from initiateParameters and are thus fixed.
% The variables of the objective function are z=(x,p).
% objFct((x,p)) = 0.5||x*-x||^2
% Constraints are of the form x_{n+1} - approx(x_n,p) = 0
% and nonnegativity constraints for p


%%%%%
% Set up fixed values of the system
%%%%%

% actual parameters, corresponding to the actual state
p_opt = pFromParams(params);
fix_params = nonOptParams(params);
start_params = pToParams(p_0,fix_params);

%%
% ODE is to be solved on the intervall [0,T]
T = 14;
t_span = [0,T];

%%
% initial value of [s,theta,sd,thetad]
x_0 = [1; -pi*0.7/2; 0; 0];

%%
% torques defined as piecewise constant functions
% Later, in Discretization, use supporting points
u1 = [  -15000;
        -15000;
        -25000;
        -25000;
        -30000;
        -30000;
        -30000;
        -30000;
        -30000;
        -25000;
        -30000;
        -30000;
        -20000;
        -20000;
        -20000;
        -20000;
        -20000;
        -20000;
        -20000;
        -20000;
        -15000;
        -15000;
        -10000;
        -15000;
        -15000;
        -15000;
        -15000;
        -20000;
        -20000;
        -20000];
u2 = [  -15000;
        -20000;
        -20000;
        -10000;
        -10000;
        -10000;
        -10000;
        -10000;
        -15000;
        -15000;
        -20000;
        -15000;
        -10000;
        -12000;
        -12000;
        -12000;
        -12000;
        -12000;
        -6000;
        -7000;
        -8000;
        -9000;
        -15000;
        -20000;
        -25000;
        -30000;
        -30000;
        -30000];


tau_B1 = piecwConst(u1,T);
tau_B2 = piecwConst(u2,T);



%%%%%
% Set up discretization of the system
%%%%%

%%
% time points for discretization:
% t_pts = [0=t_0,t_1,...,t_N=T] = [0,T/N,2T/N,...,(N-1)T/N,T]
% mesh = [t_1-t_0,t_2-t_1,...,t_N-t_{N-1}] = [T/N,...,T/N]
t_pts = (0:N)*T/N;       % time points of discretization
mesh = ones(1,N)*T/N;    % time steps of discretization

%%
% Solve the ODE
% evaluate the solution on the given discretization
% x_opt is then the discrete state for the given parameters
% x_opt will be used for the optimization as reference state
% x_opt ∊ R^4x(N+1)
% x_opt_vect ∊ R^4*(N+1)
sol_ode = solveLagrODE(params,x_0,t_span,tau_B1,tau_B2);
x_opt = deval(sol_ode,t_pts);
x_opt_vect = x_opt(:);


%%
% discrete torque values to avoid repetitive evaluation
tau_B1_val = tau_B1(t_pts);
tau_B2_val = tau_B2(t_pts);
u = [tau_B1_val(:)';
     tau_B2_val(:)'];


% Number of x values, x=[s,theta,sd,thetad]
Nx = 4*(N+1);
Np = 10;

    function [f,g] = objFctFct(z)
        % Tracking type objective function, with derivative if requested
        % z = (x,p)
        % x = [s0,θ0,sd0,θd0,s1,θ1,sd1,θd1,...]
        % f     = 1/2 |x_opt - x|^2
        % f_dx  = -(x_opt - x)
        % f_dp  = 0

        if nargout < 2
            f = 0.5*norm(x_opt_vect - z(1:Nx),2)^2;
        end

        if nargout >= 2
            f = 0.5*norm(x_opt_vect - z(1:Nx),2)^2;
            g = zeros(Nx+Np,1);
            g(1:Nx) = -(x_opt_vect - z(1:Nx));
        end
    end

    %{
    function [f,g] = objFctFct(z)
        % Tracking type objective function, with derivative if requested
        % z = (x,p)
        % x = [s0,θ0,sd0,θd0,s1,θ1,sd1,θd1,...]
        % only compare s and θ
        % f     = 1/2 |x_opt - x|^2
        % f_dx  = -(x_opt - x)
        % f_dp  = 0

        if nargout < 2
            f = 0.5*norm(x_opt_vect(1:4:end) - z(1:4:Nx),2)^2 + ...
                0.5*norm(x_opt_vect(2:4:end) - z(2:4:Nx),2)^2;

        end

        if nargout >= 2
            f = 0.5*norm(x_opt_vect(1:4:end) - z(1:4:Nx),2)^2 + ...
                0.5*norm(x_opt_vect(2:4:end) - z(2:4:Nx),2)^2;
            g = zeros(Nx+Np,1);
            g(1:4:Nx) = -(x_opt_vect(1:4:end) - z(1:4:Nx));
            g(2:4:Nx) = -(x_opt_vect(2:4:end) - z(2:4:Nx));
        end
    end
    %}



    function [c,ceq,GC,GCeq] = constrFctFct(z)
        % Constraint function for multiple shooting using explicit Euler
        % z = (x,p)
        % x_n = [s_n,θ_n,sd_n,θd_n]
        %
        % c  <= 0 :
        % -p                           <= 0
        %
        % ceq = 0 :
        % x_0 - x_0_given               = 0
        % x_{n+1} - x_n - h*rhs(x_n,p)  = 0
        currparams = pToParams(z(Nx+1:end),params);

        % inequality constraint for parameters, c<=0
        c = -z(1+Nx:end);
        %GC = zeros(Nx+Np,Np);
        GC = sparse(Nx+Np,Np);
        GC(1+Nx:end,:) = -eye(Np);

        % equality constraint from multiple shooting
        ceq = zeros(Nx,1);
        % the initial point is also used for the first values of x
        ceq(1:4) = z(1:4) - x_0(:);

        if nargout <= 2
            % no derivative information
            for i=0:N-1
                % x_{n+1} - x_n - h*rhs(x_n,p)  = 0
                % x_0     has indizes 1, ..., 4
                % x_n     has indizes 1+4*i    , ..., 4*(i+1)
                % x_{n+1} has indizes 1+4*(i+1), ..., 4*(i+2)
                rhs = rhsODE(z(1+4*i:4*(i+1)),u(:,i+1),currparams);
                ceq(1+4*(i+1):4*(i+2)) = z(1+4*(i+1):4*(i+2)) - z(1+4*i:4*(i+1)) - mesh(i+1)*rhs;
            end
        end

        if nargout >=3
            % provide derivative information
            %GCeq = zeros(Nx+Np,Nx);
            GCeq = sparse(Nx+Np,Nx);
            GCeq(1:4,1:4) = eye(4);
            for i=0:N-1
                % x_{n+1} - x_n - h*rhs(x_n,p)  = 0
                % x_0     has indizes 1, ..., 4
                % x_n     has indizes 1+4*i    , ..., 4*(i+1)
                % x_{n+1} has indizes 1+4*(i+1), ..., 4*(i+2)
                [rhs,rhs_dp,rhs_dx] = rhsODE(z(1+4*i:4*(i+1)),u(:,i+1),currparams);
                ceq(1+4*(i+1):4*(i+2)) = z(1+4*(i+1):4*(i+2)) - z(1+4*i:4*(i+1)) - mesh(i+1)*rhs;

                % deriving x_{n+1} part
                GCeq(1+4*(i+1):4*(i+2),1+4*(i+1):4*(i+2)) = eye(4);

                % deriving -x_n -h*rhs part w.r.t. x
                GCeq(1+4*i:4*(i+1),1+4*(i+1):4*(i+2)) = -eye(4) -mesh(i+1)*rhs_dx';

                % deriving -x_n -h*rhs part w.r.t. p
                GCeq(Nx+1:end,1+4*(i+1):4*(i+2)) = -mesh(i+1)*rhs_dp';
            end
        end
    end

%{
%%
% Approximation function
% calculate state approximation depending on the parameters
% choose from explicit Euler, RK order 2, 3 and 4 with index 1/2/3/4
approxFct = chooseApprox(fix_params,x_opt,u,mesh,rkOrder);

%%
% Objective function of tracking type
% Compare approximation from p to optimal approximation
% 
% |x_opt - x(p_opt)| gives an approximation error.
% To avoid this error in the objective function, one may
% use |x(p_opt) - x(p)| as the difference of the approximations
if ref == 1
    x_ref = x_opt;
else
    x_ref = approxFct(p_opt);
end
objFct = trackingTypeFct(x_ref, approxFct);
%}

objFct = @(z) objFctFct(z);
constrFct = @(z) constrFctFct(z);

% initial starting point
% use provided p_0, calculate x by forward Euler solution of ODE
x = forwEuler(start_params,x_0,u,mesh);
x_start_vect = x(:);
z_0 = [x_start_vect;p_0];

end
