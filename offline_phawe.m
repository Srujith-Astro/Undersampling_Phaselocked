clc; clear; close all;

%% =========================
% 1. LOAD DATA
%% =========================
data = readmatrix('25hz time.txt');
t = data(:,1);        % time [s]
s = data(:,2);        % signal [g]
dt = mean(diff(t));
fs = 1/dt;

%% =========================
% 2. REMOVE DC + BANDPASS FILTER
%% =========================
s = s - mean(s);

f0 = 14;   % Hz (known excitation)
bw = 2;    % bandwidth ± Hz
[b,a] = butter(3, [(f0-bw)/(fs/2), (f0+bw)/(fs/2)], 'bandpass');
s_filt = filtfilt(b, a, s);

%% =========================
% 3. PHASE ESTIMATION (Hilbert)
%% =========================
z      = hilbert(s_filt);
phi    = unwrap(angle(z));       % continuous unwrapped phase
phi_mod = mod(phi, 2*pi);        % wrapped to [0, 2π]

%% =========================
% 4. TARGET PHASE
%% =========================
phi_target = 0;    % radians — change to any value in [0, 2π]

%% =========================
% 5. FIND PHASE CROSSINGS (robust wrap-around detection)
%% =========================
% Shift phase so the target becomes the wrap point, then detect the jump
phi_shifted = mod(phi - phi_target, 2*pi);

% A crossing occurs where the shifted phase wraps (large negative jump)
cross_idx = find(diff(phi_shifted) < -pi);

% Interpolate exact crossing time for sub-sample accuracy
t_trigger = zeros(length(cross_idx), 1);
for i = 1:length(cross_idx)
    k  = cross_idx(i);
    t1 = t(k);   t2 = t(k+1);
    p1 = phi_shifted(k);
    p2 = phi_shifted(k+1) + 2*pi;   % unwrap the jump for interpolation
    t_trigger(i) = t1 + (0 - p1) * (t2 - t1) / (p2 - p1);
end

%% 6. REMOVE SPURIOUS TRIGGERS
fprintf('Triggers before cleaning: %d\n', length(t_trigger));

if length(t_trigger) > 2
    dt_trig   = diff(t_trigger);
    T_expected = 1/f0;                          % use known frequency, not mean
    valid      = true(size(t_trigger));
    valid(2:end) = abs(dt_trig - T_expected) < 0.3 * T_expected;  % ±30% of 1/14s
    t_trigger  = t_trigger(valid);
else
    warning('Not enough trigger points.');
end

fprintf('Triggers after cleaning:  %d\n', length(t_trigger));  % should be ~420
%% =========================
% 7. INTERPOLATE PHASE & AMPLITUDE AT TRIGGER TIMES
%% =========================
phi_at_trigger = interp1(t, phi_mod,  t_trigger);   % phase [0, 2π]
amp_at_trigger = interp1(t, s_filt,   t_trigger);   % filtered signal amplitude [g]

%% =========================
% 8. EXPORT CSV
%% =========================
T = table(t_trigger, phi_at_trigger, amp_at_trigger, ...
    'VariableNames', {'TriggerTime_s', 'Phase_rad', 'Amplitude_g'});

writetable(T, 'trigger_times.csv');
disp('CSV saved: trigger_times.csv');
disp(head(T));   % preview first few rows in Command Window

%% =========================
% 9. PLOT — VERIFY RESULTS
%% =========================
figure('Name', 'Phase Trigger Verification', 'NumberTitle', 'off');

subplot(3,1,1)
plot(t, s_filt, 'b', 'LineWidth', 0.8); hold on
plot(t_trigger, amp_at_trigger, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6)
xlabel('Time (s)'); ylabel('Amplitude (g)')
title('Filtered signal with trigger points')
legend('Signal', 'Triggers'); grid on

subplot(3,1,2)
plot(t, phi_mod, 'Color', [0.2 0.6 0.2], 'LineWidth', 0.8); hold on
yline(phi_target, 'r--', 'LineWidth', 1.5)
plot(t_trigger, phi_at_trigger, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6)
xlabel('Time (s)'); ylabel('Phase (rad)')
title('Wrapped phase [0, 2π] with crossings')
legend('Phase', 'Target', 'Triggers'); grid on

subplot(3,1,3)
if length(t_trigger) > 1
    dt_trig = diff(t_trigger);
    plot(t_trigger(2:end), dt_trig, 'k.-')
    yline(1/f0, 'r--')
    xlabel('Time (s)'); ylabel('\Deltat trigger (s)')
    title(sprintf('Inter-trigger interval (expected = %.4f s)', 1/f0))
    grid on
end

sgtitle('Phase Trigger Verification — 14 Hz', 'FontSize', 13)