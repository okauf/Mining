function [A_inv_b_dx] = lagrddDx(params)
% Create a function for calculating the derivative of the solution of the lagrangean 2nd order equation
% Input:
%   params      fixed parameters of the system
% Output:
%   A_inv_b_dx  function for calculating d(A_inv*b)/dx, depending on x and u
%
% The lagrangean formalism results in a system of equations of the type
%       A*[sdd; thetadd] = b
%       A*y = b
%
% The aim is to calculate dy/dx, x = [x1,x2,x3,x4]
% It holds:
%       d(A*y)/dxi = dA/dxi*y + A*dy/dxi = db/dxi
% Thus:
%       A*dy/dxi = db/dxi - dA/dxi*y,
% and dy/dxi can be calculated as the solution of a linear equation.


function y_dx = A_inv_b_dx_fct(x,u)
    % Calculate the solution of A*dy/dx = db/dx - dA/dx*y

    % second derivatives of rope 1 which are explicitely needed here
    % they depend on x, so also contribute with a derivative term
    [a_s,b_s,c_s,a_theta,b_theta,c_theta] = derderRopeLength(params, x(1), x(2), x(3), x(4)); 
    [a_s_dx,b_s_dx,c_s_dx,a_theta_dx,b_theta_dx,c_theta_dx] = derderRopeLengthDx(params, x(1), x(2), x(3), x(4)); 

    % generalized forces
    % they depend on x, so also contribute with a derivative term
    [Q_s,Q_theta] = genForces(params, x(1), x(2), x(3), x(4), u(1), u(2));
    [Q_s_dx,Q_theta_dx] = genForcesDx(params, x(1), x(2), x(3), x(4), u(1), u(2));

    % kinetic and potential energy derivatives
    % they depend on x, so also contribute with a derivative term
    [ dTds, dVds, dTdtheta, dVdtheta ] = kinPotE( params, x(1), x(2), x(3), x(4));
    [ dTds_dx, dVds_dx, dTdtheta_dx, dVdtheta_dx ] = kinPotEDx( params, x(1), x(2), x(3), x(4));


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
    % derivative of A w.r.t. every component in x
    matA_dx = zeros(2,2,4);

    % x1
    matA_dx(:,:,1) = [ a_s_dx(1)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2)),  b_s_dx(1)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2));
                       a_theta_dx(1)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2)), params.M2*(4*params.l5 - 2*params.l4 + 4*x(1)) + params.M1*(2*params.l5 + 2*x(1)) + b_theta_dx(1)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2))];

    % x2
    matA_dx(:,:,2) = [ a_s_dx(2)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2)), b_s_dx(2)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2));
                       a_theta_dx(2)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2)), b_theta_dx(2)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2))];

    % x3
    matA_dx(:,:,3) = [ a_s_dx(3)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2)), b_s_dx(3)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2));
                       a_theta_dx(3)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2)), b_theta_dx(3)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2))];

    % x4
    matA_dx(:,:,4) = [ a_s_dx(4)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2)), b_s_dx(4)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2));
                       a_theta_dx(4)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2)), b_theta_dx(4)*(params.I_B1/(2*params.r_B1^2) + params.I_P1/(2*params.r_P1^2))];

    %%
    % derivative of b w.r.t. every component in x
    % b has mostly the same structure,
    % with one individual term for x1, x3 and x4 each
    bVec_dx = [ Q_s_dx(:)' + dTds_dx(:)' - dVds_dx(:)' - ...
                    0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*c_s_dx(:)';
                Q_theta_dx(:)' + dTdtheta_dx(:)' - dVdtheta_dx(:)' - ...
                    0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*c_theta_dx(:)' ];

    % x1 addition
    bVec_dx(:,1) = bVec_dx(:,1) + [ 0;
                                    -2*x(3)*x(4)*(params.M1 + 2*params.M2)];

    % x3 addition
    bVec_dx(:,3) = bVec_dx(:,3) + [ 0;
                                    -x(4)*(4*params.M2*(params.l5 - params.l4/2 + x(1)) + 2*params.M1*(params.l5 + x(1)))];

    % x4 addition
    bVec_dx(:,4) = bVec_dx(:,4) + [ 0;
-x(3)*(4*params.M2*(params.l5 - params.l4/2 + x(1)) + 2*params.M1*(params.l5 + x(1)))];

    %%
    % Calculate dy/dxi as solution of
    %   A*dy/dxi = db/dxi - dA/dxi*y
    % Maybe only one big linear equation? %%%%%
    y_dx = zeros(2,4);
    for i=1:4
        y_dx(:,i) = linsolve( matA, bVec_dx(:,i) - matA_dx(:,:,i)*y);
    end


end

% nested functions cannot be returned, thus define an anonymous function
A_inv_b_dx = @(x,u) A_inv_b_dx_fct(x,u);


end
