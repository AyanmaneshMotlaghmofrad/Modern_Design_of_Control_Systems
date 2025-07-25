% lab 09 problem 01
clc
clear variables
close all
format compact

%% system specifications
s = tf('s')
K_min = 9;
K_max = 16;
Kn = (K_max+K_min)/2;

p1_min = 0.55;
p1_max = 1.05;
p1_n = (p1_min + p1_max)/2;

p2_min = 1.9;
p2_max = 3.1;
p2_n = (p2_min + p2_max)/2;

Ga = 0.095;
Gs = 1;
Gf = 1;
Gpn = Kn/(s*(1+ s/p1_n)*(1+ s/p2_n));
%% Weigthing functions
WS_inv = 0.15*s*(1+s/0.02)/(1+(2*0.707*s/0.426)+(s/0.426)^2);
omega = logspace(-4,4,5000);

figure(1)
bodemag(omega,WS_inv,'b')
xlabel('frequency')
grid on
WT_inv = 1.05/(1+ (2*0.658*s/2.8) + (s/2.8)^2);



%% Wu and WT
Wu = 1.1*((s + 0.3)*(s + 1.7))/((s+ 0.64)*(s + 3.1));

figure(2)
bodemag(omega,Wu)
hold on
bodemag(omega,1/WT_inv,'r')
xlabel('frequency')
grid on


% max (Wu,WT) = WT = W2

%% W1 and W2
W1 = 1/WS_inv;

lambda = 0.01;

W1mod = minreal(zpk(W1*(s/(s+lambda))),1e-4)

figure(3)
bodemag(omega,WS_inv)
hold on, grid on
bodemag(omega,1/W1mod)
xlabel('frequency')

% W2
W2mod = tf(1/1.05);

%% Designing the controller

[Am,Bm,Cm,Dm] = linmod('gen_plant');
M = ltisys(Am,Bm,Cm,Dm); % generalized plant
M = sderiv(M,2,[1/2.8 1]);
M = sderiv(M,2,[1/2.8 1]);

[gopt,Gcmod1] = hinflmi(M,[1 1], 0, 1e-2, [0 0 0]);

%% 

[Ac,Bc,Cc,Dc] = ltiss(Gcmod1);
sys = ss(Ac,Bc,Cc,Dc);
Gcmod = minreal(zpk(sys),1e-4);

Gc = minreal(Gcmodz,1e-4)

%% Loop functions

Ln = Gcmod*Ga*Gpn*Gs*Gf;
omega_Ln = logspace(-6,6,10000);
figure
nichols(Ln,omega_Ln)











