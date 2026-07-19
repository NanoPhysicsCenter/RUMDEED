% Kristinn Torfason
% 06.03.15
% Plot the charge density

clear all
close all

q_e = -1.6E-19;

r_tip = 20E-10;
eta_1 = 0.9;
theta_tip = acos(eta_1);
a_foci = r_tip/(sin(theta_tip)*tan(theta_tip));

sigma_max = 3.523*q_e/a_foci^2;

% Read the data
data = importdata('run_50x50_sumint/charge_graph.dt');

% Find the number of grid points
N = max(data(:, 1)); % x
M = max(data(:, 2)); % y

% Read the data into vectors
x = data(:, 3);
y = data(:, 4);
z = data(:, 7);

% Create the matrices for the grid and results
X = zeros(M, N);
Y = zeros(M, N);
Z = zeros(M, N);

% Convert it to the format Matlab wants
for i = 1:N
    for j = 1:M
        X(i, j) = x(i + (j - 1)*N);
        Y(i, j) = y(i + (j - 1)*N);
        if (abs(z(i + (j - 1)*N)) > 0.00125)
            Z(i, j) = 0.0;
        else
            Z(i, j) = z(i + (j - 1)*N);
        end
    end
end

% Plot the data
surf(X, Y, Z);
%axis vis3d
grid on
axis([min(x) max(x) min(y) max(y) 0 1.5E-3]);
view(0, 0);
xlabel('x [nm]');
ylabel('y [nm]');
zlabel('z [nm]');