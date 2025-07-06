clc
clear variables
close all
format compact

%% Control system setting

s =tf('s');
Gp = 25/(s^3 + 3.3*s^2 + 2*s);
Gs = 1;
Ga = 0.095;
Kd = 1;
Gf = 1/(Gs*Kd);

%% After requirement analysis;
Kc = 6; 

nu = 0;


Tp_max = 1.05;
Sp_max = 1.36;
wc_min = 0.66;
wc_max = 1.4;

%% Initial Loop function

Lin= Kc/(s^nu)*Gp*Ga*Gf*Gs; 

figure,
myngridst(Tp_max,Sp_max)
nichols(Lin)
title("Lin")
%% selection of wc desired
wcd =1;

%% positive pole
p = 0.2;
% Rp = (1+(s/10*p))/(1-(s/p));
Rp = (1 + s/100);
L= Lin * Rp;
figure,
myngridst(Tp_max,Sp_max);
nichols(L);
title("modified")

[z, p, k] = zpkdata(L/(1+L),'v')