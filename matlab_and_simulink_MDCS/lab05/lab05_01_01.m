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
wc_min = 0.66;
wc_max = 1.4;
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
bodemag(S2,{0.0001,100},'--b'),
bodemag(S_at_0,'--g')
bodemag(tf(Sp_max),'--r')
bodemag(tf(10^(-32/20)),{0.0001,0.02},'k')
grid on
% weight function

W_S = (S_star0*s)*(1+s/0.003)*(1+s/0.73)/((1+s/0.03)*(1+1.6*s/0.79+(s/0.79)^2))
bodemag(W_S,{0.0001,100},'-k')

%% T weight function
T2 = wn2^2/(s^2 + 2*zeta2*wn2*s + wn2^2);

figure,
bodemag(T2,{0.0001,100},'--b'),
hold on, 

bodemag(tf(Tp_max),'--r')
bodemag(tf(10^(-46/20)),{40,100},'k')
grid on

% butterworth weighting function for T
W_T =Tp_max/((1+ (2*0.71*s/2.7) + (s/2.7)^2))
 
bodemag(W_T,{0.0001,100},'-k')
