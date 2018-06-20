% KT
% Planar ion
% 15.05.2018

clear all
close all

global q V_0 d epsilon_0

%--------------------------------------------------------------------------
% Physical constants
q = +1.602176565E-19; % Electron charge
epsilon_0 = 8.854187817E-12; % \epsilon_0

% System parameters
w_theta = 4.7;
V_0 = 500.0;
d = 1000.0E-9;
R = 250.0E-9;

% Ion
x_a = 0.0;
y_a = 0.0;
z_a = linspace(0.5E-9, 1000.0E-9, 10000);


% Planar integral
R_plane = linspace(0.0, R, 10000);
Delta_R = R_plane(2) - R_plane(1);

for i = 1:length(z_a)
    p_d(i) = z_a(i);
    E = E_tot(R_plane, 0.0, p_d(i));
    FN_st = FN_J(w_theta, -E);
    I(i) = 2*pi*sum(R_plane.*FN_st)*Delta_R;
end

I_vac = Planar_I(V_0, w_theta, R, d);

figure()
semilogy(p_d/1E-9, I, [p_d(1)/1E-9, p_d(end)/1E-9], [I_vac, I_vac]);
xlabel('d [nm]')
ylabel('log(I) [log(A)]')
xlim([0.5, 12.5])
%ylim([10^-45, 10^6])

%figure()
%plot(R_plane/1E-9, E);


function [E] = E_tot(x, y, z_b)
  global q V_0 d epsilon_0
  
  z = 0.0; % in the plane
  E_vac = -V_0/d;
  
  fac = -q./(2*pi*epsilon_0);

  tmp_dis_a = ( (x - 0.0).^2 + (y - 0.0).^2 + (z - z_b).^2);
  
  E_ic = fac*z_b./sqrt(tmp_dis_a).^3;
  E = E_ic + E_vac;
  
end