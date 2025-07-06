% MATLAB code for magnitude and phase response of 1 / (1 + jω/p)

% Define the normalized frequency vector (logarithmic scale)
w_over_p = logspace(-2, 2, 500);  % Frequency range from 10^-2 to 10^2

% Calculate the transfer function response
H = 1 ./ (1 + 1i * w_over_p);

% Magnitude in dB
magnitude_db = 20 * log10(abs(H));

% Phase in degrees (for positive real pole, phase increases)
phase_deg = angle(H) * (180 / pi);  % Convert radians to degrees
phase_deg = -phase_deg;  % Invert the phase for a positive real pole

% Create the figure for plotting
figure;

% Plot the magnitude response
subplot(2, 1, 1);
semilogx(w_over_p, magnitude_db, 'k', 'LineWidth', 2);
title('Frequency Response of 1/(1 + jw/p) for a Positive Real Pole');
ylabel('Magnitude (dB)');
grid on;

% Plot the phase response
subplot(2, 1, 2);
semilogx(w_over_p, phase_deg, 'k', 'LineWidth', 2);
ylabel('Phase (degrees)');
xlabel('Normalized Frequency \omega/p');
grid on;
