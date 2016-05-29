function [test_dlr1d2ds,test_dlr1d2dtheta] = Check_DerRopeLength( params ) %,test_Der_dlr1d2dsd,test_Der_dlr1d2dthetad -> Der wrt t

syms s theta sd thetad sdd thetadd

[ dlr1d2ds, dlr1d2dtheta, Der_dlr1d2dsd, Der_dlr1d2dthetad ] = DerRopeLength( params, s, theta, sd, thetad, sdd, thetadd );

[ grads, gradt, gradsd, gradtd] = SymComRope( params );

test_dlr1d2ds = isequaln(dlr1d2ds,grads);
test_dlr1d2dtheta = isequaln(dlr1d2dtheta,gradt);

% test_Der_dlr1d2dsd = isequaln(Der_dlr1d2dsd,gradsd);
% test_Der_dlr1d2dthetad = isequaln(Der_dlr1d2dthetad,gradtd);

end

