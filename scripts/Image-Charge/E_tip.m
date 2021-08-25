function [E] = E_tip(xi, V_0, a, eta_1)
%E_TIP Summary of this function goes here
%   Detailed explanation goes here

eta_2 = 0.0;

fac_sqrt = 1.0 ./ ( sqrt(xi.^2 - eta_1^2)*sqrt(1 - eta_1^2 ) );
fac_log = 1.0 / log( (1+eta_1)/(1-eta_1) * (1-eta_2)/(1+eta_2) );

E = abs( 2*V_0/a* fac_sqrt * fac_log );
end

