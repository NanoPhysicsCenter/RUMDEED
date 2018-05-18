function [I] = FN_I(V_0, R_r, h_r, d, w_theta)
%FN_I Summary of this function goes here
%   Detailed explanation goes here

max_xi = h_r/d + 1;
a = sqrt(d^2*R_r^2/(h_r^2+2*d*h_r) + d^2);
eta_1 = - d / a;

max_x = a*sqrt(max_xi^2-1)*sqrt(1-eta_1^2)*1;
x_tip = linspace(-max_x, max_x, 10001);
xi = sqrt(x_tip.^2/(a^2*(1-eta_1^2)*1) + 1);

%xi = linspace(1.0, max_xi, 10001*2);

F = E_tip(xi, V_0, a, eta_1);

FN_st = FN_J(w_theta, F);

Delta_xi = abs(xi(2) - xi(1));

I = 2*pi*a^2*sqrt(1-eta_1^2)*sum(FN_st.*sqrt(xi.^2-eta_1^2))*Delta_xi/2;
end

