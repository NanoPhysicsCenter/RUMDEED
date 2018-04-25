
close all
clear variables

eta_p = linspace(-0.971, -0.975900072948533, 10000);


for i = 1:length(eta_p)
    %[FN(i), I(i), p_d(i)] = Sphere_IC_2D_v3_with_field_function(eta_p(i));
    [FN(i), I(i), p_d(i), E(i)] = Sphere_IC_2D_v4_function(eta_p(i));
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
plot(p_d/1E-9, E)
xlabel('d [nm]')
ylabel('E [V/m]')