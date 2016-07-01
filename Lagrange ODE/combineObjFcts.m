function newObjFct = combineObjFcts(objFcts)
% Combine several functions with value and gradients output

    function [ z, z_dp ] = combined(p)
        if nargout < 2
            % only output function value
            z = objFcts{1}(p);
            for i=2:length(objFcts)
                z = z + objFcts{i}(p);
            end
        end

        if nargout >= 2
            % output function value and gradient
            [z,z_dp] = objFcts{1}(p);

            for i=2:length(objFcts)
                [zi,z_dpi] = objFcts{i}(p);
                z = z + zi;
                z_dp = z_dp + z_dpi;
            end
        end
    end

newObjFct = @(p) combined(p);

end
