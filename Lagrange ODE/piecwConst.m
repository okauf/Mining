function fct = piecwConst(y,T)
% Create a piecewise constant function from input vector.
% Input:
%   y       vector which contains function values for every interval
%   T       endpoint of domain [0,T] of the function
% Output:
%   fct     piecewise constant function @(t)
%
% fct: [0,T] -> R is piecewise constant on equally long intervalls.
% Let n = length(y), then
% fct(t) = y(1) for t ∊ [0,T/n),
% fct(t) = y(2) for t ∊ [T/n,2T/n),
% ...
% fct(t) = y(n) for t ∊ [(n-1)T/n,T].

n = length(y);

    function z = pwFct(t)
        i = floor(t*n/T);
        i = min(i,n-1); % in case t == T
        z = y(i+1);
    end

fct = @(t) pwFct(t);
end
