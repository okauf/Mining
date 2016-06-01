function [ alpha, l, pos_load ] = ExcavatorConfig( lr1, lr2, params )
% given: rope length of pulleys
% computation of configuration of excavator

l = lr2 - params.l1;

% alpha = acos((params.l2^2+(rope1-params.l3)^2-(l+params.l5)^2)/2*params.l2*(rope1-params.l3));
alpha = acos(-((lr2-params.l1+params.l5)^2-params.l2^2-(lr1-params.l3)^2)/(2*params.l2*(lr1-params.l3)));

MP = [cos(params.ang_base);sin(params.ang_base)]*params.l3;

if abs(alpha + params.ang_base - pi/2) > 1e-6

    if params.ang_base + alpha > pi/2
        phi = alpha-pi/2+params.ang_base;
        pos_load = MP + [sin(phi);-cos(phi)]*(lr1-params.l3); %% rechts
    else
        phi = pi/2-alpha-params.ang_base;
        pos_load = MP - [sin(phi);cos(phi)]*(lr1-params.l3); %%links
    end

else
    pos_load = MP + [0;-lr1+params.l3]; %% senkrecht runter
end

end

