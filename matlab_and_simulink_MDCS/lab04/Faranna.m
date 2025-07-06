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
%% S weigth function
close all
wn2 = 0.91
zeta2 = 0.59

S2 = s*(s + 2*zeta2*wn2)/(s^2 + 2*zeta2*wn2*s + wn2^2) ;
S_at_0 = S_star0*s^(nu + p); %system type 1 

% plotting
figure(1)
hold on, 
bodemag(S2,{0.0001,100},'--b'), %second order
bodemag(S_at_0,'--g')
bodemag(tf(Sp_max),'--r')
bodemag(tf(10^(-32/20)),{0.0001,0.02},'k')
grid on
% weight function
z1 = 0.0025; %the first zero
p1 = 0.022; %the first pole
z2 = 0.6; %the second zero
p2 = 0.77;
W_S_negative = S_star0*s^(nu + p)*(1+s/z1)*(1+s/z2)/((1+s/p1)*(1+(1.42*s/p2)+(s/p2)^2))
bodemag(W_S_negative,{0.0001,100},'-k')

%% T weight function
T2 = wn2^2/(s^2 + 2*zeta2*wn2*s + wn2^2);

bodemag(T2,{0.0001,100},'--b'),


bodemag(tf(Tp_max),'--r')
bodemag(tf(10^(-46/20)),{40,100},'k')

% butterworth weighting function for T

p3 = 2.762% the first two poles
W_T_negative =Tp_max/(1+(1.42*s/p3)+(s/p3)^2);
bodemag(W_T_negative,{0.0001,100},'-k')
