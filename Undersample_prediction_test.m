%% Initialization
clc; clear; close all;

% Parameters
f = 5;                  % Signal frequency (Hz)
fs = 2000;              % Sampling rate (Reduced for simulation clarity)
duration = 2;           % Total simulation time
t_total = 0:1/fs:duration;
omega = 2*pi*f;

% Simulation Signal (with a random initial phase)
true_phase_rad = rand * 2 * pi; 
raw_signal = sin(omega * t_total + true_phase_rad);

% Window Settings
window_size = (1/f) * 2 * fs; % 2 cycles window
buffer = zeros(1, floor(window_size));
t_ref = (0:length(buffer)-1)/fs;

% Reference Signals (Unit Normalized)
sin_ref = sin(omega * t_ref);
cos_ref = cos(omega * t_ref);

%% Real-Time Simulation Loop
figure('Color', 'w');
h_plot = animatedline('LineWidth', 2);
axis([0 duration -1.2 1.2]);
grid on; title('Real-Time Phase Tracking');

for i = 1:length(raw_signal)
    % 1. Update Buffer (Shift left and add new sample)
    buffer = [buffer(2:end), raw_signal(i)];
    
    % 2. Pre-processing (Remove DC offset)
    proc_buffer = buffer - mean(buffer);
    
    % 3. Correlation (Dot Product)
    a = dot(proc_buffer, sin_ref); % In-phase component
    b = dot(proc_buffer, cos_ref); % Quadrature component
    
    % 4. Phase Estimation
    % atan2 returns [-pi, pi]. We convert to degrees.
    est_phase_deg = rad2deg(atan2(b, a));
    
    % 5. Visualization (Every 50 samples to save CPU)
    if mod(i, 50) == 0
        addpoints(h_plot, t_total(i), raw_signal(i));
        fprintf('Current Time: %.2fs | Estimated Phase: %.2f°\n', t_total(i), est_phase_deg);
        drawnow limitrate;
    end
end