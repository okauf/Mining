function [A_inv_b_dp] = lagrddDp(params)
% Create a function for calculating the derivative of the solution of the lagrangean 2nd order equation
% Input:
%   params      fixed parameters of the system
% Output:
%   A_inv_b_dp  function for calculating d(A_inv*b)/dp, depending on x and u
%
% The lagrangean formalism results in a system of equations of the type
%       A*[sdd; thetadd] = b
%       A*y = b
%
% The aim is to calculate dy/dp, p = [p1,...,pk,...,p10]
% It holds:
%       d(A*y)/dpk = dA/dpk*y + A*dy/dpk = db/dpk
% Thus:
%       A*dy/dpk = db/dpk - dA/dpk*y,
% and dy/dpk can be calculated as the solution of a linear equation.


function y_dp = A_inv_b_dp_fct(x,u)
    % Calculate the solution of A*dy/dp = db/dpk - dA/dp*y

    % second derivatives of rope 1 which are explicitely needed here
    % they do not depend on p
    [a_s,b_s,c_s,a_theta,b_theta,c_theta] = derderRopeLength(params, x(1), x(2), x(3), x(4)); 

    % generalized forces
    % they depend on p, so also contribute with a derivative term
    [Q_s,Q_theta] = genForces(params, x(1), x(2), x(3), x(4), u(1), u(2));
    [Q_s_dp,Q_theta_dp] = genForcesDp(params, x(1), x(2), x(3), x(4), u(1), u(2));

    % kinetic and potential energy derivatives
    % they depend on p, so also contribute with a derivative term
    [ dTds, dVds, dTdtheta, dVdtheta ] = kinPotE( params, x(1), x(2), x(3), x(4));
    [ dTds_dp, dVds_dp, dTdtheta_dp, dVdtheta_dp ] = kinPotEDp( params, x(1), x(2), x(3), x(4));


    % Calculate matrix A of the equation
    % those formulas originate from Der_dTdsd and Der_dTdthetad
    matA = [params.M1 + params.M2 + params.I_B2/params.r_B2^2 + params.I_P2/params.r_P2^2 + 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*a_s, ...
            0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*b_s;
            0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*a_theta, ...
            params.M1*(x(1) + params.l5)^2 + params.M2*(0.25*params.r_M2^2 + 1/12*params.l4^2 + 2*(x(1) + params.l5 - params.l4/2)^2) + 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*b_theta];


    % Calculate vector b of the equation
    % explicit formulas originate from Der_dTdsd and Der_dTdthetad
    bVec = [Q_s + dTds - dVds ...
                - 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*c_s;
            Q_theta + dTdtheta - dVdtheta ...
                - (2*params.M1*(x(1) + params.l5) + 4*params.M2*( x(1) + params.l5 - 0.5*params.l4))*x(3)*x(4) ...
                - 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*c_theta];

    y = linsolve( matA, bVec );


    %%
    % derivative of A w.r.t. every parameter in p
    matA_dp = zeros(2,2,10);

    % M1
    matA_dp(:,:,1) = [ 1,   0;
                        0, (x(1) + params.l5)^2];
    % M2
    matA_dp(:,:,2) = [ 1,   0;
                       0,   0.25*params.r_M2^2 + 1/12*params.l4^2 + 2*(x(1) + params.l5 - 0.5*params.l4)^2];
    % I_B1
    matA_dp(:,:,3) = [ a_s/(2*params.r_B1^2),      b_s/(2*params.r_B1^2);
                       a_theta/(2*params.r_B1^2),  b_theta/(2*params.r_B1^2)];
    % I_B2
    matA_dp(:,:,4) = [ 1/params.r_B2^2,    0;
                        0,          0];
    % I_P1
    matA_dp(:,:,5) = [ a_s/(2*params.r_P1^2),      b_s/(2*params.r_P1^2);
                       a_theta/(2*params.r_P1^2),  b_theta/(2*params.r_P1^2)];
    % I_P2
    matA_dp(:,:,6) = [ 1/params.r_P2^2,    0;
                        0,                 0];
    % mu_B1, mu_B2, mu_P1, mu_P2 contribute with zero matrices,
    % already in matA_dp by default


    %%
    % derivative of b w.r.t. every parameter in p
    % b has mostly the same structure,
    % only M1, M2, I_B1 and I_P1 have additional terms
    bVec_dp = [Q_s_dp(:)' + dTds_dp(:)' - dVds_dp(:)';
                    Q_theta_dp(:)' + dTdtheta_dp(:)' - dVds_dp(:)'];

    % M1 addition
    bVec_dp(:,1) = bVec_dp(:,1) + [ 0;
                                    -2*x(3)*x(4)*(x(1) + params.l5)];
    % M2 addition
    bVec_dp(:,2) = bVec_dp(:,2) + [ 0;
                                    -4*x(3)*x(4)*(x(1) + params.l5 - 0.5*params.l4)];
    % I_B1 addition
    bVec_dp(:,3) = bVec_dp(:,3) + [ -c_s/(2*params.r_B1^2);
                                    -c_theta/(2*params.r_B1^2)];
    % I_P1 addition
    bVec_dp(:,4) = bVec_dp(:,4) + [ -c_s/(2*params.r_P1^2);
                                    -c_theta/(2*params.r_P1^2)];

    %%
    % Calculate dy/dpk as solution of
    %   A*dy/dpk = db/dpk - dA/dpk*y
    % Maybe only one big linear equation? %%%%%
    y_dp = zeros(2,10);
    for k=1:10
        y_dp(:,k) = linsolve( matA, bVec_dp(:,k) - matA_dp(:,:,k)*y);
    end


end

% nested functions cannot be returned, thus define an anonymous function
A_inv_b_dp = @(x,u) A_inv_b_dp_fct(x,u);


end
