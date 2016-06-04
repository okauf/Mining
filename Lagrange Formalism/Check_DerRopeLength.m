
% Comparison of implemented derivatives with derivatives calculated using
% the symbolic toolbox (apart from derivatives wrt time)
 
syms s(t) theta(t) sd(t) thetad(t) sdd(t) thetadd(t)

[ dlr1d2ds, dlr1d2dtheta, Der_dlr1d2dsd, Der_dlr1d2dthetad ] = DerRopeLength( params, s, theta, sd, thetad, sdd, thetadd );

[ grads, gradt, Der_gradsd, Der_gradtd] = SymComRope( params );

%Check successful if test variables equal 0
test_dlr1d2ds = simplify(dlr1d2ds-grads);
test_dlr1d2dtheta = simplify(dlr1d2dtheta-gradt);

test_Der_dlr1d2dsd = simplify(Der_dlr1d2dsd-Der_gradsd);
test_Der_dlr1d2dthetad = simplify(Der_dlr1d2dthetad-Der_gradtd);


