%clear all
close all

V_0 = 2E3;
d = 1000E-9;
L = linspace(10, 1000, 10) * 1E-9;
h = 0.0;
x = 0.15;

q_e = -1.602176565E-19;
m_e = 9.10938291E-31;
epsilon_0 = 8.854187817E-12;
h_bar = 1.054571726E-34;

a_FN = q_e^2/(16*pi^2*h_bar);
b_FN = -4/(3*h_bar) * sqrt(-2*m_e*q_e);
l_const = -q_e / (4*pi*epsilon_0);
w_theta = 2.0;

E_vac = V_0 / d ;
F = 1;

N = 50;
%J_keep = 1:N;
%F_keep = 1:N;

J_L = 1:length(L);
F_L = 1:length(L);
LD_C = 1:length(L);

J_CL = 4/9*epsilon_0*sqrt(-2*q_e/m_e)*V_0^(3/2)/d^2;

for k = 1:length(L)
    LD_C(k) = L(k)/(2*d);
    %DL_C = d/L(k);
    
    J_old = 0.0;
    J_keep(1:N) = 0.0;
    F_keep(1:N) = 0.0;
    for i = 1:N
        F = E_vac * F;
        l = l_const * F / w_theta^2;
        if (l > 1)
            %disp('Warning l > 1');
            disp('L');
            disp(L(k)/1E-9);
            disp('F');
            disp(F);
            error('Error l > 1')
        end
        t_y = 1 + l*(1/9 - 1/18*log(l));
        v_y = 1 - l + 1./6 * l * log(l);
        %t_y = 1;
        %v_y = 1;
        
        elec_supply = a_FN * F^2 / (w_theta * t_y^2);
        esc_prob = exp(b_FN * w_theta^(3/2) * v_y / F);
        
        J_new = elec_supply * esc_prob / J_CL;
        J = x*J_new + (1-x)*J_old;
        
        fun = @(y, z) sqrt(z) ./ ((y.^2 + z.^2).*sqrt(LD_C(k).^2 + y.^2 + z.^2));
        E_z = quad2d(fun, -LD_C(k), LD_C(k), 0.0, 1, 'AbsTol', 1E-3, 'RelTol', 1E-3, 'Singular', true, 'MaxFunEvals', 100000);
        E_z = J/(9*pi) * 2*LD_C(k) * E_z;

%         fun = @(y, z) sqrt(z) ./ ((y.^2 + DL_C^2.*z.^2).*sqrt(0.25 + y.^2 + DL_C^2.*z.^2));
%         E_z = quad2d(fun, -0.5, 0.5, 0.0, 1, 'AbsTol', 1E-3, 'RelTol', 1E-3, 'Singular', true, 'MaxFunEvals', 100000);
%         E_z = J*DL_C/(9*pi) * E_z;
        
        F = 1 - 2*E_z;
        %E_z = E_z * V_0/d;
        %F = E_vac + 2*E_z;
        %if (F > 0.0)
        %    F = 0.0;
        %    disp('Warning F > 0');
        %else
        %    F = -F;
        %end
        
        F_keep(i) = F;
        J_keep(i) = J;
        J_old = J;
    end
    
    if abs((J_keep(N) - J_keep(N-1))) > 1E-3
        disp('Warning error > 1E-3');
        disp('L');
        disp(L(k)/1E-9);
        disp('Error');
        disp(abs((J_keep(N) - J_keep(N-1))));
    end
    
    J_L(k) = J;
    F_L(k) = F;
end

figure();
plot(L/1E-9, J_L*J_CL);
title('Value');
xlabel('L [nm]');
ylabel('J [A/m^2]');

figure();
plot(LD_C, J_L);
title('Scaled CL');
xlabel('L/2d');
ylabel('J/J_{CL}^{1D}');

figure();
plot(L/1E-9, F_L);
title('Surface field');
xlabel('L [nm]');
ylabel('F');

% figure();
% plot(LD_C, J_L*J_CL/7.4543E12);
% title('Scaled SC');
% xlabel('L/2d');
% ylabel('J/J_{SC}^{1D}');

%fileID = fopen('output.dt', 'w');
%fprintf(fileID, '%12.8f %12.8f\n', [L/1E-9; J_L]);
%fclose(fileID);