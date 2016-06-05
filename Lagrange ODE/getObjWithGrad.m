function [objWithGrad] = getObjWithGrad(objFct,objGrad)

    function [f,g] = objWithGradFct(p)
        % Calculate objective value f
        f = objFct(p);
        
        if nargout > 1 % gradient required
            g = objGrad(p);
        end
    end

objWithGrad = @(p) objWithGradFct(p);

end