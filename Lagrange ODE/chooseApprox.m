function [ approxFct ] = chooseApprox(params, x_opt, u, mesh, i)
% Select one of the approximation functions
% Input:
%   params      fix parameters, not containing optimizable parameters
%   x_opt       optimal reference value, 4xN
%   u           torque values on supporting points of both functions
%   mesh        time steps of discretization
%   i           index of approx function
%               1:  explicit Euler
%               2:  RK of order 2
%               3:  RK of order 3
%               4:  RK of order 4
% Output:
%   approxFct   function depending on parameters to calculate approximation

if i == 1
    approxFct = explEulerApprox(params, x_opt, u, mesh);
elseif i == 2
    approxFct = rk2Approx(params, x_opt, u, mesh);
elseif i == 3
    approxFct = rk3Approx(params, x_opt, u, mesh);
else
    approxFct = rk4Approx(params, x_opt, u, mesh);
end
