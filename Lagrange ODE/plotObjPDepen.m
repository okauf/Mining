function plotObjPDepen(p_opt, objFct, d)
% Plot dependency of the objective function w.r.t. to each parameter
% Input:
%   p_opt   optimal parameters
%   objFct  objective function
%   d       relative deviation of the parameter values
%
% Fix all parameters except one. Plot the objective function value
% for different values of this one free parameter.
% The values for this free parameter are in
%   [(1-d)s,(1+d)s],
% where d is the deviation and s the optimum value of the parameter.

% names of the parameters
p_name = {'$M_1$', '$M_2$', '$I_{B_1}$', '$I_{B_2}$', ...
          '$I_{P_1}$', '$I_{P_2}$', '$\mu_{B_1}$', ...
          '$\mu_{B_2}$', '$\mu_{P_1}$', '$\mu_{P_2}$'};

block_name = {'Mass', 'Inertia', 'Friction'};

    function [ y ] = obj_i(p_i,i)
        p = p_opt;
        p(i) = p_i;
        y = objFct(p);
    end

    function plotParamBlock(m,n,ka,kb,b)
        % m,n   dimensions of the subplot
        % ka,kb index range for parameters to plot
        % b     block number

        figure;
        for i=ka:kb
            obj = @(p_i) obj_i(p_i,i);
            p_i = linspace((1-d)*p_opt(i),(1+d)*p_opt(i),10);
            y = zeros(length(p_i),1);
            for j=1:length(p_i)
                y(j) = obj(p_i(j));
            end
            subplot(m,n,i-ka+1);
            plot(p_i,y);
            title(p_name(i),'interpreter','latex');
        end
    end

plotParamBlock(1,2,1,2,1);
plotParamBlock(2,2,3,6,2);
plotParamBlock(2,2,7,10,3);


%{
figure;
for i=1:length(p_opt)
    obj = @(p_i) obj_i(p_i,i);
    p_i = linspace((1-d)*p_opt(i),(1+d)*p_opt(i),10);
    y = zeros(length(p_i),1);
    for j=1:length(p_i)
        y(j) = obj(p_i(j));
    end
    subplot(5,2,i);
    plot(p_i,y);
end
%}

end
