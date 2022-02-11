clear all
close all

%run('/home/kristor/physconst.m');

r_tip = 20E-10;
eta_1 = 0.9;
theta_tip = acos(eta_1);
a_foci = r_tip/(sin(theta_tip)*tan(theta_tip));
%a_foci= 94.7E-10;

length_scale = 1E-9;

N = 100;
M = 100;

x_start = -100.0*length_scale;
x_end   = -1.0*x_start;
y_start = -100.0*length_scale;
y_end   = -1.0*y_start;

x_l = x_end - x_start;
y_l = y_end - y_start;

x_inc = x_l / (N-1);
y_inc = y_l / (M-1);

for i = 1:N
  x = x_start + (i-1)*x_inc;
  
  for j = 1:M
    y = y_start + (j-1)*y_inc;
    
    phi_p = atan2(y, x);
    xi_p = 1.0/a_foci * 1.0/sqrt(1-eta_1^2) * sqrt(x^2 + y^2 + a_foci^2*(1-eta_1^2));
  end

end

x = -5E-9;
y = 2E-9;
xi_p = 1.0/a_foci * 1.0/sqrt(1-eta_1^2) * sqrt(x^2 + y^2 + a_foci^2*(1-eta_1^2))

xi = linspace(1.0, 2.0, 100);
phi = linspace(0, 2*pi, 100);

x = a_foci * sqrt(xi.^2 - 1) * sqrt(1-eta_1^2) .* cos(phi);
y = a_foci * sqrt(xi.^2 - 1) * sqrt(1-eta_1^2) .* cos(phi);
z = a_foci * xi * eta_1;
