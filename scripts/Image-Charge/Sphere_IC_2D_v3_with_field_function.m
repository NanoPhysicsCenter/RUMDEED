function [FN, I, p_d, E_vac, E_tot, E_ic] = Sphere_IC_2D_v3_with_field_function(xi_p, eta_p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Kristinn Torfason
% 24.04.2015
% Sphere approximation of image charge effect V3
% One sphere at location nearest to the charge
%
% Here we find the point on the tip closest to the electron. There we
% calculate the radius of curvature of the tip. We then create a sphere
% with that radius, that passes through this point and is orthogonal to a
% line from the electron to the point on that tip. This sphere is then
% use to calculate the image charge partner. See chapter 2.2 in Jackson.

%--------------------------------------------------------------------------
% Physical constants
q = +1.602176565E-19; % Electron charge
m_e = 9.10938291E-31; % Electron mass
epsilon_0 = 8.854187817E-12; % \epsilon_0
h_bar = 1.054571726E-34;
w_theta = 4.7; % eV

%--------------------------------------------------------------------------
% Constants the define the tip geometry

%Base radius
R = 250.0E-9; % [nm]

% Height
h = 500E-9; % [nm]

% Gap spacing, measured from the top of the tip.
d = 1000.0E-9; % [nm]
d_plane = d + h;

% Voltage
V_0 = 500.0;

max_xi = h/d + 1;
a = sqrt(d^2*R^2/(h^2+2*d*h) + d^2);
eta_1 = - d / a;
theta = acos(d/a);
r_tip = a*sin(theta)*tan(theta);
shift_z = abs(a*eta_1*max_xi);

%disp('a')
%disp(a/1.0E-9)

%
max_x = a*sqrt(max_xi^2-1)*sqrt(1-eta_1^2)*1;
x_tip = linspace(-max_x, max_x, 10001);
%x_tip = 2.5719e-09;
%x = 0.0;
y_tip = 0.0;

xi = sqrt(x_tip.^2/(a^2*(1-eta_1^2)*1) + 1);
z_tip = a .* xi .* eta_1;

%--------------------------------------------------------------------------
% Particle
% Location in prolate spheroidal coordinates
% https://en.wikipedia.org/wiki/Prolate_spheroidal_coordinates#Alternative_definition
% (sigma, tau, phi) in Wiki, here (xi, eta, phi).

% Height of the particle, ξ [1, inf].
% xi = 1 is the top of the tip.
%xi_a1 = 1.020;
%xi_a1 = 1.0;
xi_a1 = xi_p;

% Distance from the tip, η [0, -1].
% Should not be lower than eta_1. That would be inside the tip.
% Anything above 0 would be above the absorption plane.
%eta_a1 = -0.97515; % 1 nm, if xi = 1.02.
%eta_a1 = -0.9745; % 2 nm, if x = 1.02.
%eta_a1 = -0.9745;
eta_a1 = eta_p;

% Rotation of the particle around the tip.
phi_a1 = 0.0;

% (x, y, z) coordinate of the particle.
x_a = a * sqrt(xi_a1^2 - 1) * sqrt(1 - eta_a1^2) * cos(phi_a1);
y_a = a * sqrt(xi_a1^2 - 1) * sqrt(1 - eta_a1^2) * sin(phi_a1);
z_a = a * xi_a1 * eta_a1;

%disp('Particle x [nm]')
%disp(x_a1/1.0E-9)
%disp('Particle y [nm]')
%disp(y_a1/1.0E-9)
%disp('Particle z [nm]')
%disp((z_a1+shift_z)/1.0E-9)


[x_0, y_0, z_0] = Find_closest_point(x_a, y_a, z_a, eta_1, a);

% Particle distance from tip
p_d = sqrt((x_0 - x_a)^2 + (y_0 - y_a)^2 + (z_0 - z_a)^2);
%disp('Particle distance from tip [nm]')
%disp(p_d / 1.0E-9)

xi_0 = abs(z_0/(a*eta_1));
R_sphere = abs(a/eta_1 * ( (xi_0^2 - eta_1^2).^(3/2) ) / (sqrt(1 - eta_1^2)));

n(1) = eta_1/sqrt(1-eta_1^2)*x_0/sqrt(x_0^2 + y_0^2 + a^2*(1-eta_1^2));
n(2) = eta_1/sqrt(1-eta_1^2)*y_0/sqrt(x_0^2 + y_0^2 + a^2*(1-eta_1^2));
n(3) = -1;

n = n/norm(n);

x_c = n(1).*R_sphere + x_0;
y_c = n(2).*R_sphere + y_0;
z_c = n(3).*R_sphere + z_0;

% Distance from center of sphere (x_c, y_c, z_c) to particle (x_a, y_a, z_a).
a_b = sqrt((x_a - x_c).^2 + (y_a - y_c).^2 + (z_a - z_c).^2);

% Distance of image charge from center of sphere.
b = r_tip^2/a_b;

% Vector from center of sphere to particle. Unit length.
n_a(1, :) = (x_a - x_c) ./ a_b;
n_a(2, :) = (y_a - y_c) ./ a_b;
n_a(3, :) = (z_a - z_c) ./ a_b;

% Coordinates of image charge partner.
x_b = x_c + n_a(1, :).*b;
y_b = y_c + n_a(2, :).*b;
z_b = z_c + n_a(3, :).*b;

%--------------------------------------------------------------------------
% Calculate the electric field.

tmp_dis_a = ( (x_tip - x_a).^2 + (y_tip - y_a).^2 + (z_tip - z_a).^2 );
tmp_dis_b = ( (x_tip - x_b).^2 + (y_tip - y_b).^2 + (z_tip - z_b).^2 );

% Size of the charge of the image charge partner.
q_ic = -r_tip/a_b*q;

% x
E_x1 = q/(4*pi*epsilon_0)    * (x_tip - x_a)./sqrt(tmp_dis_a).^(3);
E_x2 = q_ic/(4*pi*epsilon_0) * (x_tip - x_b)./sqrt(tmp_dis_b).^(3);
E_x = E_x1 + E_x2;

% y
E_y1 = q/(4*pi*epsilon_0)    * (y_tip - y_a)./sqrt(tmp_dis_a).^(3);
E_y2 = q_ic/(4*pi*epsilon_0) * (y_tip - y_b)./sqrt(tmp_dis_b).^(3);
E_y = E_y1 + E_y2;

% z
E_z1 = q/(4*pi*epsilon_0)    * (z_tip - z_a)./sqrt(tmp_dis_a).^(3);
E_z2 = q_ic/(4*pi*epsilon_0) * (z_tip - z_b)./sqrt(tmp_dis_b).^(3);
E_z = E_z1 + E_z2;

% Magnitude of the field.
E_ic = sqrt(E_x.^2 + E_y.^2 + E_z.^2);

[Ex_vac, Ey_vac, Ez_vac] = Tip_Field(x_tip, y_tip, z_tip, eta_1, a, V_0);

E_x = E_x + Ex_vac;
E_y = E_y + Ey_vac;
E_z = E_z + Ez_vac;

E_new = sqrt(E_x.^2 + E_y.^2 + E_z.^2);
E_vac = sqrt(Ex_vac.^2 + Ey_vac.^2  + Ez_vac.^2);

E_tot = E_new;

%x_s = x_tip/1E-9;
%sigma_s = E_new*epsilon_0/1E-6;
%save('data_sphere_z-1475-xi-1_0.mat', 'x_s', 'sigma_s', '-v7')

a_FN = q^2/(16.0*pi^2*h_bar); % A eV V^{-2}
b_FN = -4.0/(3.0*h_bar) * sqrt(2.0*m_e*q); % eV^{-3/2} V m^{-1}
l_const = q / (4.0*pi*epsilon_0); % eV^{2} V^{-1} m

l = l_const * E_new / w_theta^2;
l(l>1.0) = 1.0;
l(l<0.0) = 0.0;
v_y = 1 - l + 1/6*l.*log(l);
t_y = 1 + l.*(1/9 - 1/18*log(l));

FN_st = a_FN ./ (t_y.^2 .* w_theta) .* E_new.^2 .* exp(b_FN .* w_theta^(3/2) .* v_y ./ E_new);

%FN = FN_st(length(FN_st)/2);
[FN, m_l] = max(FN_st);

Delta_xi = abs(xi(2) - xi(1));
I = 2*pi*a^2*sqrt(1-eta_1^2)*sum(FN_st.*sqrt(xi.^2-eta_1^2))*Delta_xi/2;
end

