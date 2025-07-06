%% homework 05
% Modern Design of Control Systems
clc
clear variables
close all
format compact

%% defining the system

s = tf('s');

Gp = 25/(s^3+3.3*s^2+2*s);
Ga = 0.095;
Gs = 1;
Kd = 1; %requirement 1
Gf = 1/(Kd*Gs);

Tp_max = 1.05;
Sp_max = 1.36;
wc_min = 0.66;%non serve
wc_max = 1.4;%non serve
S_star0 = 0.15;
nu = 0;
p = 1;
M_LF = -32;

%% weighting function on S W_S_inv

close all

wn2 = 0.91;
zeta2 = 0.59;

S2 = s*(s + 2*zeta2*wn2)/(s^2 + 2*zeta2*wn2*s + wn2^2);
S_0 = s^(nu+p)*S_star0;

figure(1)
grid on
bodemag(S2,{0.0001,1000},'--b')
hold on
bodemag(S_0,{0.0001,10},'--g')
bodemag(tf(Sp_max),{0.0001,1000},'--r')
bodemag(tf(10^(M_LF/20)),{0.0001,0.02})

z1 = 0.003;
p1 = 0.03;
z2 = 0.8;
p2 = sqrt(z1*z2*Sp_max/(S_star0*p1));
zeta =0.85;
W_S_inv = S_0 * (1 + s/z1)*(1 + s/z2)/((1+(2*zeta*s/p2)+(s/p2)^2)*(1+s/p1));
bodemag(W_S_inv,{0.0001,1000},'-k')

%% Weighting function on T: W_T
close 2
T2 = wn2^2/(s^2 + 2*zeta2*wn2*s + wn2^2);
p3 = 6.18;

figure(2)
bodemag(T2,{0.001,1000},'--b')
hold on, grid on
bodemag(tf(Tp_max),{0.001,1000},'--r')
bodemag(tf(10^(-32/20)),{40,1000},'--k')

W_T_inv = Tp_max/(1+(2*0.707*s/p3)+(s/p3)^2);
bodemag(W_T_inv,{0.001,1000},'k')

%%






