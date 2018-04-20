% Kristinn Torfason
% 24.04.2015
% Sphere approximation of image charge effect V4
% Sphere location fixed at tip

close all
clear variables

q = +1.602176565E-19; % Electron charge
epsilon_0 = 8.854187817E-12; % \epsilon_0

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

max_xi = h_r/d + 1;
a = sqrt(d^2*R_r^2/(h_r^2+2*d*h_r) + d^2);
eta_1 = - d / a;

theta = acos(d/a);
r_tip = a*sin(theta)*tan(theta);
eta_2 = 0.0; % Defines the plane of absorption
eta = eta_1; % Defines the tip
shift_z = abs(a*eta*max_xi);

%--------------------------------------------------------------------------
% Particle
% Position in prolate spheroidal coordinates
xi_a = 1.020;
eta_a = -0.8745;
phi_a = 0.0;

% Position of the particle in x,y,z coordinates
x_a = a * sqrt(xi_a^2 - 1) * sqrt(1 - eta_a^2) * cos(phi_a);
y_a = a * sqrt(xi_a^2 - 1) * sqrt(1 - eta_a^2) * sin(phi_a);
z_a = a * xi_a * eta_a;


%-------------------------------------------------------------------------
% Caclulate x, y and z values for the tip
% Max x value, should be the same as R_r (base radius)
max_x = a*sqrt(max_xi^2-1)*sqrt(1-eta_1^2)*1;
x_tip = linspace(-max_x, max_x, 1000);
%x_tip = 2.5719e-09;
%x = 0.0;
y_tip = 0.0;

xi = sqrt(x_tip.^2/(a^2*(1-eta_1^2)*1) + 1);
z_tip = a .* xi .* eta_1;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw the tip, sphere and charge.
figure(1);
hold on
axis equal
plot(x_tip/1E-9, z_tip/1E-9); % Tip

t = linspace(0, 2*pi, 100);
x_sp = x_c + g_fac*r_tip*cos(t);
z_sp = z_c + g_fac*r_tip*sin(t);

plot(x_sp/1E-9, z_sp/1E-9); % Sphere

plot(x_a/1E-9, z_a/1E-9, 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 2); % Charge
plot(x_b/1E-9, z_b/1E-9, 'ok', 'MarkerFaceColor', 'r', 'MarkerSize', 2); % Partner
plot([-max_x/1E-9, max_x/1E-9], [0, 0]); % Top plane
xlabel('x [nm]')
ylabel('z [nm]')

%--------------------------------------------------------------------------
% Calculate the electric field.

tmp_dis_a = ( (x_tip - x_a).^2 + (y_tip - y_a).^2 + (z_tip - z_a).^2 );
tmp_dis_b = ( (x_tip - x_b).^2 + (y_tip - y_b).^2 + (z_tip - z_b).^2 );

% Size of the charge of the image charge partner.
q_ic = -r_tip/a_b*q;

% x
E_x1 = q/(4*pi*epsilon_0)    * (x_a - x_tip)./sqrt(tmp_dis_a).^(3);
E_x2 = q_ic/(4*pi*epsilon_0) *  (x_b - x_tip)./sqrt(tmp_dis_b).^(3);
E_x = E_x1 + E_x2;

% y
E_y1 = q/(4*pi*epsilon_0)    * (y_a - y_tip)./sqrt(tmp_dis_a).^(3);
E_y2 = q_ic/(4*pi*epsilon_0) * (y_b - y_tip)./sqrt(tmp_dis_b).^(3);
E_y = E_y1 + E_y2;

% z
E_z1 = q/(4*pi*epsilon_0)   * (z_a - z_tip)./sqrt(tmp_dis_a).^(3);
E_z2 = q_ic/(4*pi*epsilon_0) * (z_b - z_tip)./sqrt(tmp_dis_b).^(3);
E_z = E_z1 + E_z2;

% Magnitude of the field.
E_new = sqrt(E_x.^2 + E_y.^2 + E_z.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw the field as a function of position on the tip.
figure(2);
plot(x_tip/1E-9, E_new);
xlabel('x [nm]');
ylabel('E [V/m]');

%x_nmv4 = x_tip / 1E-9;
%sigma_muv4_xi1_030 = E_new*epsilon_0/1E-6;

%%save('V4.mat', 'x_nmv4', 'sigma_muv4_xi1_000')
%save('V4.mat', 'sigma_muv4_xi1_030', '-append')