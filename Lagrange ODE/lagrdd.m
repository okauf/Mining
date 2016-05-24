function [A_inv_b] = lagrdd(params)
% Create a function for calculating the solution of the lagrangean 2nd order equation
% Input:
%   params      fixed parameters of the system
% Output:
%   A_inv_b     function for calculating A_inv*b, depending on x and u
%
% The lagrangean formalism results in a system of equations of the type
%       A*[sdd; thetadd] = b
% A depends on the state and parameters
% b depends on the state, control and parameters
% Thus, A_inv_b has input x,u


function y = A_inv_b_fct(x,u)
    % Calculate the solution of A*y = b

    % second derivatives of rope 1 which are explicitely needed here
    [a_s,b_s,c_s,a_theta,b_theta,c_theta] = derderRopeLength(params, x(1), x(2), x(3), x(4)); 

    % generalized forces
    [Q_s,Q_theta] = genForces(params, x(1), x(2), x(3), x(4), u(1), u(2));

    % kinetic and potential energy derivatives
    [ dTds, dVds, dTdtheta, dVdtheta ] = kinPotE( params, x(1), x(2), x(3), x(4));


    % Calculate matrix A of the equation
    % those formulas originate from Der_dTdsd and Der_dTdthetad
    matA = [params.M1 + params.M2 + params.I_B2/params.r_B2^2 + params.I_P2/params.r_P2^2 + 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*a_s, ...
            0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*b_s;
            0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*a_theta, ...
            params.M1*(x(1) + params.l5)^2 + params.M2*(0.25*params.r_M2^2 + 1/12*params.l4^2 + 2*(x(1) + params.l5 - params.l4/2))^2 + 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)];


    % Calculate vector b of the equation
    % explicit formulas originate from Der_dTdsd and Der_dTdthetad
    bVec = [Q_s + dTds - dVds ...
                - 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*c_s;
            Q_theta + dTdtheta - dVdtheta ...
                - (2*params.M1*(x(1) + params.l5)^2 + 4*params.M2*( x(1) + params.l5 - 0.5*params.l4))*x(3)*x(4) ...
                - 0.5*(params.I_B1/params.r_B1^2 + params.I_P1/params.r_P1^2)*c_theta];


    y = linsolve( matA, bVec );

end

% nested functions cannot be returned, thus define an anonymous function
A_inv_b = @(x,u) A_inv_b_fct(x,u);


end
