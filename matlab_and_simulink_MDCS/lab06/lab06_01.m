clc 
clear variables
close all
format compact

%% defining the system

s = tf('s');
Kmin = 9;
Kmax = 16;
Kn = (Kmin + Kmax)/2;

p1min = 0.55;
p1max = 1.05;
p1n = (p1min+p1max)/2;

p2min = 1.9;
p2max = 3.1;
p2n = (p2min+p2max)/2;

Gpn = Kn/(s*(1+s/p1n)*(1+s/p2n));
Gs = 1;
Ga = 0.095;

%% Plotting the class of TF
close all

omega = logspace(-2,2,5000);
Gpn_mag = squeeze(freqresp(Gpn, omega));

mag_max = -Inf*ones(length(omega),1); % Initialization for storing max magnitude

figure(1);
hold on;
grid on;

for K = linspace(Kmin, Kmax, 10)
    for p1 = linspace(p1min, p1max, 10)
        for p2 = linspace(p2min, p2max, 10)
            Gp = K / (s * (1 + s/p1) * (1 + s/p2));
            Gp_mag = squeeze(freqresp(Gp, omega));
            delta_m_mag = abs(Gp_mag ./ Gpn_mag - 1);
            loglog(omega, delta_m_mag); % Using loglog directly
            
            % getting the max
            mag_max = max(mag_max,delta_m_mag);
        end
    end
end

plot(omega, mag_max, 'LineWidth', 2, 'Color', 'r');
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title('Frequency Response Comparison');

%% Fitting Wu
figure(2)
magg = vpck(mag_max,omega);
Wu = fitmag(magg);
[A,B,C,D] = unpck(Wu);

[Z_w,P_w,K_w] = ss2zp(A,B,C,D);
%%
sys = ss(A,B,C,D);
Wu= zpk(tf(sys))
figure(1)
Wu_mag = abs(squeeze(freqresp(Wu, omega)));
loglog(omega,Wu_mag, 'LineWidth', 2, 'Color', 'g');
set(gca,'XScale','log')

%%


