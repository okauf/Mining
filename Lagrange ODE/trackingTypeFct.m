function trckFct = trackingTypeFct(x_opt, approxFct)
% Create a tracking type function with norm 2
% Input:
%   x_opt       optimal reference value
%   approxFct   function depending on parameters to calculate approximation
% Output:
%   trckFct     function to calculate z =  0.5|| x_opt - approxFct(p) ||


    function [ z, z_dp ] = trackingType(p)
        % Tracking type objective function, with derivative if requested
        % z     = 1/2 |x_opt - x(p)|^2
        % z_dp  = -(x_opt - x(p)'*x_dp

        if nargout < 2
            x = approxFct(p);
            % x is of size 4xN
            z = 0.5*norm(x_opt(1,:) - x(1,:),2)^2 + ...
                0.5*norm(x_opt(2,:) - x(2,:),2)^2 + ...
                0.5*norm(x_opt(3,:) - x(3,:),2)^2 + ...
                0.5*norm(x_opt(4,:) - x(4,:),2)^2;
        end

        if nargout >= 2
            [x,x_dp] = approxFct(p);
            % x    is of size 4xN
            % x_dp is of size 4xNx10
            z = 0.5*norm(x_opt(1,:) - x(1,:),2)^2 + ...
                0.5*norm(x_opt(2,:) - x(2,:),2)^2 + ...
                0.5*norm(x_opt(3,:) - x(3,:),2)^2 + ...
                0.5*norm(x_opt(4,:) - x(4,:),2)^2;

            z_dp = - (x_opt(1,:) - x(1,:))*squeeze(x_dp(1,:,:)) - ...
                     (x_opt(2,:) - x(2,:))*squeeze(x_dp(2,:,:)) - ...
                     (x_opt(3,:) - x(3,:))*squeeze(x_dp(3,:,:)) - ...
                     (x_opt(4,:) - x(4,:))*squeeze(x_dp(4,:,:));
        end
    end

trckFct = @(p) trackingType(p);

end
