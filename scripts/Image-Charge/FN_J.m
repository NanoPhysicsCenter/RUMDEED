function [J] = FN_J(w_theta, F)
%FN_J Summary of this function goes here
%   Detailed explanation goes here

%--------------------------------------------------------------------------
% Physical constants
q = +1.602176565E-19; % Electron charge
m_e = 9.10938291E-31; % Electron mass
epsilon_0 = 8.854187817E-12; % \epsilon_0
h_bar = 1.054571726E-34;

a_FN = q^2/(16.0*pi^2*h_bar); % A eV V^{-2}
b_FN = -4.0/(3.0*h_bar) * sqrt(2.0*m_e*q); % eV^{-3/2} V m^{-1}
l_const = q / (4.0*pi*epsilon_0); % eV^{2} V^{-1} m

l = l_const * F / w_theta^2;
l(l>1.0) = 1.0;
l(l<0.0) = 0.0;

v_y = 1 - l + 1/6*l.*log(l);
t_y = 1 + l.*(1/9 - 1/18*log(l));

J = a_FN ./ (t_y.^2 .* w_theta) .* F.^2 .* exp(b_FN .* w_theta^(3/2) .* v_y ./ F);
end

