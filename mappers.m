clc; clear; close all;

%% Parameters
P = 64;            % Total number of antennas
L = 16;            % Length of visibility region
theta = 0.1;       % True spatial frequency (normalized, -0.5 to 0.5)
theta_grid = linspace(-0.5, 0.5, 1000);

n_full = 0:P-1;
n_vr = 0:L-1;

%% 1. VR-aware using Dirichlet envelope formula
G_vr_formula = zeros(size(theta_grid));
for i = 1:length(theta_grid)
    delta = theta - theta_grid(i);
    num = sin(pi * L * delta);
    denom = sin(pi * delta) + 1e-12;  % avoid division by zero
    G_vr_formula(i) = abs(num / (L * denom))^2;
end

%% 2. VR-aware using steering vector correlation
a_theta = exp(1j * 2 * pi * n_vr * theta).' / sqrt(L);
G_vr_vector = zeros(size(theta_grid));
for i = 1:length(theta_grid)
    a_prime = exp(1j * 2 * pi * n_vr * theta_grid(i)).' / sqrt(L);
    G_vr_vector(i) = abs(a_prime' * a_theta)^2;
end

%% 3. Blind case: Full steering vector masked over P antennas
k = 20;  % start index of VR
mask = zeros(P, 1);
mask(k+1 : k+L) = 1;

a_theta_blind = mask .* exp(1j * 2 * pi * n_full.' * theta) / sqrt(P);
G_blind = zeros(size(theta_grid));
for i = 1:length(theta_grid)
    a_prime_blind = exp(1j * 2 * pi * n_full.' * theta_grid(i)) / sqrt(P);
    G_blind(i) = abs(a_prime_blind' * a_theta_blind)^2;
end

%% Plotting
figure;
plot(theta_grid, G_vr_formula, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_grid, G_vr_vector, 'r--', 'LineWidth', 1.5);
plot(theta_grid, G_blind, 'g-', 'LineWidth', 1.5);
legend('VR-aware (formula)', 'VR-aware (vector)', 'Blind');
xlabel('Normalized Spatial Frequency \theta''');
ylabel('Gain');
title('Beam Patterns: VR-aware vs Blind (Unnormalized)');
grid on;