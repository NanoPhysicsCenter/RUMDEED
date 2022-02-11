%
function [Ex, Ey, Ez] = Tip_Field( x_0, y_0, z_0, eta_1, a, V_0 )

eta_2 = 0.0;

xi = 1/(2.0*a) * ( sqrt(x_0.^2 + y_0.^2 + (z_0+a).^2) + sqrt(x_0.^2 + y_0.^2 + (z_0-a).^2) );
eta = eta_1;
%eta = 1/(2.0*a) * ( sqrt(x_0.^2 + y_0.^2 + (z_0+a).^2) - sqrt(x_0.^2 + y_0.^2 + (z_0-a).^2) );
theta = atan2(y_0, x_0);

pre_fac = 2.0*V_0/a * 1./(xi.^2 - eta.^2) * 1/log(( (1+eta_1)/(1-eta_1) )*( (1-eta_2)/(1+eta_2) ));
pre_fac_xy = -eta .* sqrt((xi.^2-1)./(1-eta.^2));

Ex = pre_fac.*pre_fac_xy.*cos(theta);
Ey = pre_fac.*pre_fac_xy.*sin(theta);
Ez = pre_fac.*xi;

end