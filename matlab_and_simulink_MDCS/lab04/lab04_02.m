clc
clear variables
close all
format compact

%% Defining the system and its specifications

s = tf('s');
Gp = 40/(s^2 + 3*s + 4.5);
Kp = dcgain(Gp);
p = 0;

Gs = 1;
Ga = -0.09;

%% First specification
Kd = 1;
Gf = 1/(Gs*Kd)

%% transient requirements
% S8)
s_hat = 0.08;

zeta = abs(log(s_hat))/sqrt(pi^2 + log(s_hat)^2);

Tp0 = 1/(2*zeta*sqrt(1-zeta^2))
Sp0 = (2*zeta*sqrt(2 + 4*zeta^2+2*sqrt(1+8*zeta^2)))/(sqrt(1+8*zeta^2)+4*zeta^2-1)

Tp0_dB = 20*log10(Tp0)
Sp0_dB = 20*log10(Sp0)

% S6)
tr = 2.5
wn_tr = (pi-acos(zeta))/(tr*sqrt(1-zeta^2))
wc_tr = wn_tr*sqrt(sqrt(1+4*zeta^4)-2*zeta^2)

%S7)
ts = 10;
wn_ts = (-log(0.05))/(zeta*ts)
wc_ts = wn_ts*sqrt(sqrt(1+4*zeta^4)-2*zeta^2)
