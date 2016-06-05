% Examplary calculate derivatives w.r.t. p

syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10;
p_sym = [ p1, p2, p3, p4, p5, p6, p7, p8, p9, p10];

%%
% Approximation function
approxFctSym = explEulerApproxSym(fix_params,x_opt,u,mesh);

%%
% Objective function of tracking type
objFctSym = trackingTypeFct(x_opt, approxFctSym);


%%
% Derivative of objective function w.r.t. p
% Need to uncomment a line in explEuler for this to work
% This may be very slow for N large
%y = objFct(p);
y_sym = objFctSym(p_sym);

for i=1:length(p_sym)
    yd_sym(i) = diff(y_sym,p_sym(i));
end

objGrad = @(p) double(subs(yd_sym(:),p_sym(:),p(:)));
