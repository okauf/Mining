function [ fix_params ] = nonOptParams(params)
% Evict the given parameters from optimizable parameters

fix_params = params;

% Rewrite all optimizable parameters with 0
fix_params.M1 = 0;
fix_params.M2 = 0;
fix_params.I_B1 = 0;
fix_params.I_B2 = 0;
fix_params.I_P1 = 0;
fix_params.I_P2 = 0;
fix_params.mu_B1 = 0;
fix_params.mu_B2 = 0;
fix_params.mu_P1 = 0;
fix_params.mu_P2 = 0;

end
