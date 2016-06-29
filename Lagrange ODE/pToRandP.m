function [ p_dev ] = pToRandP( p, d )
% Create a random deviation of p
% Input:
%   p           starting values
%   d           deviation ∊ (0,1)
% Output:
%   p_dev       randomly deviated values
%
% The random values are from the interval ((1-d), (1+d))*start value

% rand gives random values from (0,1)
% transformed to (-1,1)
r = 2*rand(size(p)) - 1;

% transformed to (1-d,1+d)
r = 1 + d*r;

% value p(i) is from interval (1-d,1+d)*p_ref(i)
p_dev = r.*p;


end
