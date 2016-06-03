
% Comparison of implemented derivatives with derivatives calculated using
% the symbolic toolbox

syms s theta sd thetad sdd thetadd

[ Der_dTdsd, dTds, dVds, Der_dTdthetad, dTdtheta, dVdtheta ] = KinPotE( params, s, theta, sd, thetad, sdd, thetadd );

[ gradTs,gradVs,Der_gradTsd,gradTt,gradVt,Der_gradTtd] = SymComKinPot( params );

%Check successful if test variables equal 0
test_gradTs = simplify(gradTs-dTds);
test_gradVs = simplify(gradVs-dVds);
test_gradTsd = simplify(Der_gradTsd-Der_dTdsd); %X
test_gradTt = simplify(gradTt-dTdtheta);
test_gradVt = simplify(gradVt-dVdtheta);
test_gradTtd = simplify(Der_gradTtd-Der_dTdthetad); %X


