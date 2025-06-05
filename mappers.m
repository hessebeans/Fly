clc; clear; close all;

%% Parameters
P = 64;              % Total number of antennas
L = 16;              % Size of visibility region
k = 20;              % Starting index of the visibility region
theta = 0.1;         % True spatial frequency (normalized to [-0.5, 0.5])
theta_grid = linspace(-0.5, 0.5, 1000);  % Spatial frequency axis

%% Theoretical beam gain formulas (power domain)
x = theta - theta_grid;

% Dirichlet-like envelope (sinc squared form)
env = (sin(pi * L * x) ./ (sin(pi * x) + 1e-12)).^2;

% VR-aware theoretical gain
G_vr_theory = env / L^2;

% Blind beamforming theoretical gain
G_blind_theory = env / (P * L);

%% Channel vector: same for both cases
% Construct true channel with support only in [k, k+L-1]
h = zeros(P, 1);
for n = k:(k + L - 1)
    h(n+1) = exp(1j * 2 * pi * n * theta);  % +1 for MATLAB indexing
end
h = h / norm(h);  % Normalize to unit norm

%% Correlation-based gain computation
G_vr_corr = zeros(1, length(theta_grid));
G_blind_corr = zeros(1, length(theta_grid));

for idx = 1:length(theta_grid)
    t = theta_grid(idx);

    % VR-aware beamformer (only over visible region)
    w_vr = zeros(P, 1);
    for n = k:(k + L - 1)
        w_vr(n+1) = exp(1j * 2 * pi * n * t);  % +1 for MATLAB indexing
    end
    w_vr = w_vr / norm(w_vr);  % Normalize beamformer

    % Blind beamformer (entire aperture)
    n_idx = 0:P-1;
    w_blind = exp(1j * 2 * pi * n_idx * t).' / sqrt(P);

    % Gains = |w^H * h|^2
    G_vr_corr(idx) = abs(w_vr' * h)^2;
    G_blind_corr(idx) = abs(w_blind' * h)^2;
end

%% Plot 1: Theoretical gain expressions
figure;
subplot(1, 2, 1);
plot(theta_grid, G_vr_theory, 'b-', 'LineWidth', 2); hold on;
plot(theta_grid, G_blind_theory, 'r--', 'LineWidth', 2);
legend('VR-aware (Theory)', 'Blind (Theory)', 'Location', 'Best');
xlabel('Spatial Frequency \theta''');
ylabel('Beamforming Gain');
title('Theoretical Envelope Gain');
grid on;

%% Plot 2: Gain using correlation with steering vectors
subplot(1, 2, 2);
plot(theta_grid, G_vr_corr, 'b-', 'LineWidth', 2); hold on;
plot(theta_grid, G_blind_corr, 'r--', 'LineWidth', 2);
legend('VR-aware (Corr)', 'Blind (Corr)', 'Location', 'Best');
xlabel('Spatial Frequency \theta''');
ylabel('Beamforming Gain');
title('Correlation-based Gain');
grid on;

sgtitle('Beamforming Gain Comparison: VR-aware vs Blind (1 Ray)');