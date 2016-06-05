function [ lr1, lr2 ] = stheta2lr(params, s, theta)
% Convert s and theta into lr1 and lr2


% calculate lr1-l3 via cosinus relation in general triangle
% c^2 := (lr1-l3)^2 = l2^2 + (s+l5)^2 - 2l2(s+l5)^2 * cos(phi-theta)
c_qu = params.l2^2 + (s+params.l5).^2 - 2*params.l2.*(s+params.l5).*cos(params.ang_base - theta);
c = sqrt(c_qu);
lr1 = c + params.l3;
lr2 = s + params.l1;

end