
close all
clear variables

N = 10000;

eta_p = linspace(-0.971, -0.975900072948533, N);

I_0 = zeros(1, N);
I_05 = zeros(1, N);

pd_0 = zeros(1, N);
pd_05 = zeros(1, N);

% xi_p = 1.0000;
% disp('V-3');
% [~, I_3, pd_3, E_vac3, E_tot3, E_ic3] = Sphere_IC_2D_v3_with_field_function(xi_p, -0.971);
% disp(I_3);
% disp(pd_3/1E-9);
% 
% disp('')
% disp('V-4');
% [~, I_4, pd_4, E_vac4, E_tot4, E_ic4] = Sphere_IC_2D_v4_function(xi_p, -0.971);
% disp(I_4);
% disp(pd_4/1E-9);
% 
% figure()
% hold on
% plot(E_ic3);
% plot(E_ic4);
% legend('V3', 'V4');


xi_p = 1.000;
for i = 1:length(eta_p)
    [~, I_0(i), pd_0(i), ~, ~, ~] = Sphere_IC_2D_v3_with_field_function(xi_p, eta_p(i));
    %[~, I_0(i), pd_0(i), ~, ~, ~] = Sphere_IC_2D_v4_function(xi_p, eta_p(i));
end

xi_p = 1.000;
for i = 1:length(eta_p)
    %[~, I_05(i), pd_05(i), ~, ~, ~] = Sphere_IC_2D_v3_with_field_function(xi_p, eta_p(i));
    [~, I_05(i), pd_05(i), ~, ~, ~] = Sphere_IC_2D_v4_function(xi_p, eta_p(i));
end


I_vac = FN_I(500.0, 250.0E-9, 500.0E-9, 1000.0E-9, 4.7);

figure()
semilogy(pd_0/1E-9, I_0, pd_05/1E-9, I_05, [pd_0(1)/1E-9, pd_0(end)/1E-9], [I_vac, I_vac])
xlim([0.5, 5])
xlabel('d [nm]')
ylabel('log(I) [log(A)]')
legend('V3', 'V4', 'Vac')
