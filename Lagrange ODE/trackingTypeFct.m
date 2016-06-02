function trckFct = trackingTypeFct(x_opt, approxFct)
% Create a tracking type function with norm 2
% Input:
%   x_opt       optimal reference value
%   approxFct   function depending on parameters to calculate approximation
% Output:
%   trckFct     function to calculate z =  0.5|| x_opt - approxFct(p) ||

    function z = trackingType(p)
        x = approxFct(p);
        z = 0.5*norm(x_opt(1,:) - x(1,:),2)^2 + ...
            0.5*norm(x_opt(2,:) - x(2,:),2)^2 + ...
            0.5*norm(x_opt(3,:) - x(3,:),2)^2 + ...
            0.5*norm(x_opt(4,:) - x(4,:),2)^2;

        %z = 0.5*norm(x_opt(:) - x(:),2)^2;
    end

trckFct = @(p) trackingType(p);

end
