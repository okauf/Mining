function [ Cost, Constr ] = OptimizationFunctions( x, u, mesh )
% Cost and Constraints as functions in optimization variable

Cost = @(params) CostFct(x, u, params, mesh);
Constr = @(params) ConstrFct(x, u, params, mesh);

    function [Cost] = CostFct(x, u, params, mesh)
        Cost = @(params) 0;
    end

    function[Constr] = ConstrFct(x, u, params, mesh)
        sol_ode = discretizeLagrODE( x, u, params,mesh);
        Constr = @(params) [];
    end
end

