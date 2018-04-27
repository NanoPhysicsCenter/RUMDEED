function [FN, I, p_t, E_top, v_yy, t_yy, l_yy] = Sphere_IC_2D_v4_function(eta_a)
% Kristinn Torfason
% 26.04.2018
% Sphere approximation of image charge effect V4
% Sphere location fixed at tip

%close all
%clear variables

%--------------------------------------------------------------------------
% Physical constants
q = +1.602176565E-19; % Electron charge
m_e = 9.10938291E-31; % Electron mass
epsilon_0 = 8.854187817E-12; % \epsilon_0
h_bar = 1.054571726E-34;
w_theta = 4.7; % eV

g_fac = 1.0; % Scale factor for testing (Not really used).

%--------------------------------------------------------------------------
% Parameters for the tip


%Base radius
R_r = 250.0E-9; % [nm]
% Tip height
h_r = 500.0E-9; % [nm]
% Distance from top of tip to plane
d = 1000.0E-9; % [nm]
% Distance of plane from bottom of tip
d_plane = d + h_r;
% Voltage
V_0 = 500.0;

max_xi = h_r/d + 1;
a = sqrt(d^2*R_r^2/(h_r^2+2*d*h_r) + d^2);
eta_1 = - d / a;

theta = acos(d/a);
r_tip = a*sin(theta)*tan(theta);
%eta_2 = 0.0; % Defines the plane of absorption
%eta = eta_1; % Defines the tip
%shift_z = abs(a*eta*max_xi);

%--------------------------------------------------------------------------
% Particle
% Position in prolate spheroidal coordinates
xi_a = 1.000;
%eta_a = -0.8745;
phi_a = 0.0;

% Position of the particle in x,y,z coordinates
x_a = a * sqrt(xi_a^2 - 1) * sqrt(1 - eta_a^2) * cos(phi_a);
y_a = a * sqrt(xi_a^2 - 1) * sqrt(1 - eta_a^2) * sin(phi_a);
z_a = a * xi_a * eta_a;


%-------------------------------------------------------------------------
% Caclulate x, y and z values for the tip
% Max x value, should be the same as R_r (base radius)
max_x = a*sqrt(max_xi^2-1)*sqrt(1-eta_1^2)*1;
x_tip = linspace(-max_x, max_x, 10001);
%x_tip = 2.5719e-09;
%x = 0.0;
y_tip = 0.0;

xi = sqrt(x_tip.^2/(a^2*(1-eta_1^2)*1) + 1);
z_tip = a .* xi .* eta_1;

%--------------------------------------------------------------------------
% Top of the tip
x_t = 0.0;
y_t = 0.0;
z_t = a * 1.0 * eta_1; % xi = 1, at the top

% Distance from particle to top of tip
p_t = sqrt( (x_a - x_t)^2 + (y_a - y_t)^2 + (z_a - z_t)^2 );

%--------------------------------------------------------------------------
% Center of the sphere used for the image charge approximation
x_c = 0.0;
y_c = 0.0;
z_c = a*1.0*eta_1 - r_tip*g_fac; % Located at the top of the tip


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
E_tip = sqrt(E_x.^2 + E_y.^2 + E_z.^2);

[Ex_vac, Ey_vac, Ez_vac] = Tip_Field(x_tip, y_tip, z_tip, eta_1, a, V_0);

%E_t = sqrt(E_x.^2 + E_y.^2 + E_z.^2);

E_x = E_x + Ex_vac;
E_y = E_y + Ey_vac;
E_z = E_z + Ez_vac;

E_new = sqrt(E_x.^2 + E_y.^2 + E_z.^2);
E_vac = sqrt(Ex_vac.^2 + Ey_vac.^2  + Ez_vac.^2);

N = length(E_new);

E_top = E_tip((N-1)/2+1) + E_vac((N-1)/2+1);

a_FN = q^2/(16.0*pi^2*h_bar); % A eV V^{-2}
b_FN = -4.0/(3.0*h_bar) * sqrt(2.0*m_e*q); % eV^{-3/2} V m^{-1}
l_const = q / (4.0*pi*epsilon_0); % eV^{2} V^{-1} m

%
l = l_const * E_new / w_theta^2;
l(l>1.0) = 1.0;
l(l<0.0) = 0.0;
l_yy = max(l);
v_y = 1 - l + 1/6*l.*log(l);
t_y = 1 + l.*(1/9 - 1/18*log(l));
%v_y = 1.0;
%t_y = 1.0;

v_yy = v_y((N-1)/2+1);
t_yy = t_y((N-1)/2+1);

FN_st = a_FN ./ (t_y.^2 .* w_theta) .* E_new.^2 .* exp(b_FN .* w_theta^(3/2) .* v_y ./ E_new);

%FN = FN_st(length(FN_st)/2);
%FN = max(FN_st);

% J at x = 0
%FN_st(l>1.0) = 0.0;
FN = FN_st((N-1)/2+1);

% Calculate the current
I = 2*pi*a^2*sqrt(1-eta_1^2)*sum(FN_st.*sqrt(xi.^2-eta_1^2));

end

