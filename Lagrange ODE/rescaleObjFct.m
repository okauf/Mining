function [ fct ] = rescaleObjFct( objFct, K )
% Rescale an objective function and it's gradient
% Input:
%   objFct      function returning function value and gradient
%   K           scalar to scale objFct
% Output:
%   fct         rescaled objective Function

    function [ f, g ] = rescaledObjFct(p)
        % Objective function
        % Output:
        %   f       function value
        %   g       gradient

        if nargout < 2
            f = objFct(p);
            f = K*f;
        end

        if nargout >= 2
            [f,g] = objFct(p);
            f = K*f;
            g = K*g;
        end
    end

fct = @(p) rescaledObjFct(p);

end
