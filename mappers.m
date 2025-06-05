clc; clear; close all;

%% Parameters
P = 64;                % Total number of antennas (full array)
L = 16;                % Length of visibility region (subarray)
theta_true = 0.1;      % True direction of arrival (normalized: -0.5 to 0.5)
theta_grid = linspace(-0.5, 0.5, 1000);   % Evaluation range

%% 1. VR-Aware Case
% Use subarray of length L (assumed centered at 0 for simplicity)
n_vr = 0:L-1;   % VR-aware antenna indices

% Compute beam pattern for each theta'
G_vr = zeros(size(theta_grid));
for i = 1:length(theta_grid)
    theta_prime = theta_grid(i);
    sv_target = exp(1j * 2 * pi * n_vr * theta_true).' / sqrt(L);
    sv_scan   = exp(1j * 2 * pi * n_vr * theta_prime).' / sqrt(L);
    G_vr(i) = abs(sv_scan' * sv_target)^2;
end

%% 2. Blind Case
% BS uses full array, but only L antennas are actually illuminated
% Let visible region be from index k to k+L-1
k = 20;                            % Visibility starts at antenna index 20
n_blind = 0:P-1;                   % Full array indices
mask = zeros(P,1);
mask(k:k+L-1) = 1;                % Only these antennas are active

% Channel steering vector (only over visible part)
sv_ray = exp(1j * 2 * pi * n_blind.' * theta_true) .* mask / sqrt(P);

% Compute beam pattern for each theta'
G_blind = zeros(size(theta_grid));
for i = 1:length(theta_grid)
    theta_prime = theta_grid(i);
    sv_scan = exp(1j * 2 * pi * n_blind.' * theta_prime) / sqrt(P);
    G_blind(i) = abs(sv_scan' * sv_ray)^2;
end

%% Plot 1: VR-Aware Beam Pattern
figure;
plot(theta_grid, G_vr, 'LineWidth', 2);
xlabel('Normalized spatial frequency \theta''');
ylabel('Gain |a(\theta'')^H a(\theta)|^2');
title('Beam Pattern - VR-Aware Case (L = 16)');
grid on;
ylim([0 1.1]);

%% Plot 2: Blind (VR-Unaware) Beam Pattern
figure;
plot(theta_grid, G_blind, 'LineWidth', 2);
xlabel('Normalized spatial frequency \theta''');
ylabel('Gain |a(\theta'')^H h|^2');
title('Beam Pattern - Blind Case (only L antennas active)');
grid on;
ylim([0 1.1]);