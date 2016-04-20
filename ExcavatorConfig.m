function [ alpha, l, pos_load ] = ExcavatorConfig( rope1, rope2, params )
% given: rope length of pulleys
% computation of configuration of excavator

l = rope2 - params.l1;

alpha = acos((params.l2^2+(rope1-params.l3)^2-(l+params.l5)^2)/2*params.l2*(rope1-params.l3));

%angle between excavator arms
beta = pi/2 - alpha + ...
    acos(-((rope1-params.l3-cos(alpha)*params.l2)^2-sin(alpha)^2*params.l2^2-(l+params.l5)^2)/(2*(l+params.l5)*sin(alpha)*params.l2));

pos_load = (l+params.l5)*[cos(pi-beta);sin(pi-beta)] + params.l1*[cos(alpha);sin(alpha)];

% pos_load = (rope1-params.l3)*[cos(alpha);sin(alpha)] + params.l3*[cos(params.ang_base);sin(params.ang_base)];

end

