
close all
clear variables

eta_p = linspace(-0.971, -0.975900072948533, 10000);


for i = 1:length(eta_p)
    %[FN(i), I(i), p_d(i)] = Sphere_IC_2D_v3_with_field_function(eta_p(i));
    [FN(i), I(i), p_d(i), E(i), v_y(i), t_y(i), l(i)] = Sphere_IC_2D_v4_function(eta_p(i));
end


figure()
semilogy(p_d/1E-9, FN)
%xlim([0, 5])
xlabel('d [nm]')
ylabel('J [A/m^2]')

figure()
semilogy(eta_p, I)
%xlim([0, 5])
xlabel('eta')
ylabel('I [A]')

figure()
semilogy(p_d/1E-9, I)
%xlim([0, 5])
xlabel('d [nm]')
ylabel('I [A]')

figure()
semilogy(p_d/1E-9, E)
xlabel('d [nm]')
ylabel('E [V/m]')

figure()
plot(p_d/1E-9, v_y, p_d/1E-9, t_y, p_d/1E-9, l)
xlabel('d [nm]')
ylabel('v_y/t_y/l')
legend('v_y', 't_y', 'l')

% figure()
% q = +1.602176565E-19; % Electron charge
% m_e = 9.10938291E-31; % Electron mass
% epsilon_0 = 8.854187817E-12; % \epsilon_0
% h_bar = 1.054571726E-34;
% w_theta = 4.7; % eV
% 
% a_FN = q^2/(16.0*pi^2*h_bar); % A eV V^{-2}
% b_FN = -4.0/(3.0*h_bar) * sqrt(2.0*m_e*q); % eV^{-3/2} V m^{-1}
% l_const = q / (4.0*pi*epsilon_0); % eV^{2} V^{-1} m
% 
% l = l_const * E / w_theta^2;
% v_y = 1 - l + 1/6*l.*log(l);
% t_y = 1 + l.*(1/9 - 1/18*log(l));
% %v_y = 1.0;
% %t_y = 1.0;
% 
% FN_E = a_FN ./ (t_y.^2 .* w_theta) .* E.^2 .* exp(b_FN .* w_theta^(3/2) .* v_y ./ E);
% 
% semilogy(p_d/1E-9, FN_E);
% xlabel('d [nm]')
% ylabel('J [A/m^2]')