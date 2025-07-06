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
wk = (K_max-K_min)/2;
p1_min = 0.55;
p1_max = 1.05;
p1_n = (p1_min + p1_max)/2;
wp1 = (p1_max- p1_min)/2;
p2_min = 1.9;
p2_max = 3.1;
p2_n = (p2_min + p2_max)/2;
wp2 = (p2_max - p2_min)/2
Ga = 0.095;
Gs = 1;
Gf = 1;
Gpn = Kn/(s*(1+ s/p1_n)*(1+ s/p2_n));
%% Weigthing functions
WS_inv = 0.15*s*(1+s/0.02)/(1+(2*0.707*s/0.426)+(s/0.426)^2);
omega = logspace(-4,4,5000);

WS_inv1 = 0.15*s*(1+s/0.02)/(1+(2*0.707*s/0.412)+(s/0.412)^2);

figure(1)
bodemag(omega,WS_inv,'b')
hold on
bodemag(omega,WS_inv1,'r')
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
W1 = 1/WS_inv1;

lambda = 0.01;

W1mod = minreal(zpk(W1*(s/(s+lambda))),1e-4)

figure(3)
bodemag(omega,WS_inv)
hold on, grid on
bodemag(omega,1/W1mod)
xlabel('frequency')

% W2
W2mod = tf(1/0.95);

%% Designing the controller

[Am,Bm,Cm,Dm] = linmod('gen_plant');
M = ltisys(Am,Bm,Cm,Dm); % generalized plant
M = sderiv(M,2,[1/2.8 1]);
M = sderiv(M,2,[1/2.8 1]);

[gopt,Gcmod] = hinflmi(M,[1 1], 0, 1e-2, [0 0 0]);

%% 

[Ac,Bc,Cc,Dc] = ltiss(Gcmod);
sys = ss(Ac,Bc,Cc,Dc);
Gcmod = minreal(zpk(sys),1e-4);

%% reduction of the degree of the controller
Gc = minreal(Gcmod*(s+0.009998)*(1+s/3918)/(s+1.236e-05),1e-3)

%% Loop functions

Ln = Gcmod*Ga*Gpn*Gs*Gf;
omega_Ln = logspace(-6,6,10000);
figure
nichols(Ln,omega_Ln)
hold on
Ln1 = Gc*Ga*Gs*Gf*Gpn;
nichols(Ln1,omega_Ln)

%% time response
Tn = minreal(Ln1/(1+Ln1),1e-4);
figure
step(Tn*Kd)

Sn = minreal(1/(1 + Ln1),1e-4);
dcgain(1/s*Sn)

%%
t = linspace(0,20,1000);
er = lsim(Sn,t,t);
figure
title('output error ')
plot(t,er,'LineWidth',1.5,'Color','b')

yline(0.15)


%% mu analysis
K_min = 7;
K_max = 18;
Kn = (K_max+K_min)/2;
wk = (K_max-K_min)/2;
p1_min = 0.45;
p1_max = 1.15;
p1_n = (p1_min + p1_max)/2;
wp1 = (p1_max- p1_min)/2;
p2_min = 1.5;
p2_max = 3.5;
p2_n = (p2_min + p2_max)/2;
wp2 = (p2_max - p2_min)/2;


omega = logspace(-3,3,1000);
[An,Bn,Cn,Dn] = linmod('N_scheme_RS')
N = pck(An,Bn,Cn,Dn);
Nf = frsp(N,omega);
deltaset = [-1,1;-1,1;-1,1];
mu_bounds_RS = mu(Nf,deltaset);
figure,
vplot('liv,m',mu_bounds_RS) % iv (y axis) , m (x axis)

%% performance S3 %we check at zero
Ws_mu_s3 = 1/(0.2182 * s); %s_star for requirement 3 
omega_s3 = logspace(-6,-3,1000); 

[An,Bn,Cn,Dn] = linmod('N_scheme_S3')
N = pck(An,Bn,Cn,Dn);
Nf = frsp(N,omega_s3);
deltaset = [-1,1;-1,1;-1,1;1,1];
mu_bounds_S3 = mu(Nf,deltaset);
figure,
vplot('liv,m',mu_bounds_S3) % iv (y axis) , m (x axis)


%% S4
WS_mu_S4 = tf(1/(10^(-32/20)));
omega_s4 = logspace(-6,log10(0.02),1000); 
% up to wp+
[An,Bn,Cn,Dn] = linmod('N_scheme_S4')
N = pck(An,Bn,Cn,Dn);
Nf = frsp(N,omega_s4);
deltaset = [-1,1;-1,1;-1,1;1,1];
mu_bounds_S4 = mu(Nf,deltaset);
figure,
vplot('liv,m',mu_bounds_S4) % iv (y axis) , m (x axis)

%% S5
WT_mu_s5 = tf(1/(10^(-46/20))); %1/M_T^HF
omega_s5 = logspace(log10(40),3,1000);

%
[An,Bn,Cn,Dn] = linmod('N_scheme_S5')
N = pck(An,Bn,Cn,Dn);
Nf = frsp(N,omega_s5);
deltaset = [-1,1;-1,1;-1,1;1,1];
mu_bounds_S5 = mu(Nf,deltaset);
figure,
vplot('liv,m',mu_bounds_S5) % iv (y axis) , m (x axis)
