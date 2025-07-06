clc
clear variables
close all
format compact
for i = 1:20
    figure(i)
end
    %%
K_min = 24;
K_max = 56;
K_n = (K_max + K_min)/2
Wk = (K_max - K_min)/2

b_min = 68;
b_max = 100;
b_n = (b_max + b_min)/2
Wb = (b_max - b_min)/2

c = 0.707
s = tf('s')
Gpn = zpk(K_n/(s*(1+2*c*s/b_n + (s/b_n)^2)))
kp = dcgain(s*Gpn)
Gs = 3;
Ga = -100;
Da0 = 8.5e-3;
Dp0 = 3.5e-3;
as = 1e-2;

%%
Kd = 1;
nu =1;
Gf = 1/(Kd*Gs)
M_HF = 2e-4*Gs/as
M_HF_dB = mag2db(M_HF)

%% transient
s_hat = 0.12;
zeta = -log(s_hat)/sqrt(pi^2 + log(s_hat)^2)

Tpo = 1/(2*zeta*sqrt(1-zeta^2))
Tpo_dB = mag2db(Tpo)
Spo = (2*zeta*sqrt(2+4*zeta^2+2*sqrt(1+8*zeta^2)))/(sqrt(1+8*zeta^2) + 4*zeta^2 -1)

tr = 0.007;
wn_tr = (pi-acos(zeta))/(tr*sqrt(1-zeta^2))

ts = 0.014;
wn_ts = -log(0.05)/(zeta*ts)

wn = max(wn_tr,wn_ts)

%% WS_inv
close 1
S2 = s*(s + 2*zeta*wn)/(s^2 + 2*zeta*wn*s + wn^2)

omega = logspace(-3,5,2000);


figure(1)
hold on, grid on
bodemag(omega,S2,'--k')
bodemag(omega,tf(Spo),'k')

%WS_inv for performance
p1 = 0.003;
z1 = 190;
S_star_0 =1 
p2 = sqrt(Spo*z1/(p1*S_star_0))
WS_inv = s^2*(1+s/z1)/((1+s/p1)*(1+ (2*0.8*s/p2) + (s/p2)^2))

bodemag(omega,WS_inv,'b')

% WS_inv for controller
p1 = 0.003;
z1 = 200
S_star_0 =1 
p2 = sqrt(Spo*z1/(p1*S_star_0))
WS_inv1 = s^2*(1+s/z1)/((1+s/p1)*(1+ (2*0.707*s/p2) + (s/p2)^2))

bodemag(omega,WS_inv1,'g')

%% W2
close 2
p3 = 3000*10^((Tpo_dB-M_HF_dB)/-40)

T2 = wn^2/(s^2 + 2*zeta*wn*s + wn^2);

WT_inv = Tpo/(1+ (2*0.72*s/p3)+ (s/p3)^2)

figure(2)
hold on,grid on
bodemag(omega,T2,'--k')
bodemag({3000,100000},tf(M_HF),'k')
bodemag(omega,tf(Tpo),'k')

bodemag(omega,WT_inv,'b')

%% Wu
close 4
Gpn_mag = squeeze(freqresp(Gpn,omega));
mag_max = zeros(length(omega),1);

figure(4)
hold on,grid on

for K = linspace(K_min,K_max,30)
    for b = linspace(b_min,b_max,30)
        Gp = K/(s*(1+2*c*s/b + (s/b)^2));
        Gp_mag = squeeze(freqresp(Gp,omega));
        delta = abs(Gp_mag./Gpn_mag -1);

        mag_max = max(delta,mag_max);
        plot(omega,delta)

    end
end

plot(omega, mag_max,'r','LineWidth',1.5)
set(gca,'XScale','log')
set(gca, 'YScale','log')

%%
magg = vpck(mag_max,omega);
wu = fitmag(magg);
[A,B,C,D] = unpck(wu);

Wu = zpk(tf(ss(A,B,C,D)))

%% Wu and Wt
figure(5)
hold on, grid on
bodemag(omega, Wu,'r')
bodemag(omega, 1/WT_inv,'b')
%% W1 and W2

W2 = (1+s/p3)^2/Tpo;
W2mod = tf(1/Tpo);

figure(2)
bodemag(omega,W2,'c')

W1 = 1/WS_inv1;
lambda = 5;
W1mod = minreal(W1*s^2/(s +lambda)^2,1e-3)
figure(1)
bodemag(omega,1/W1mod,'c')

%%
[Am,Bm,Cm,Dm] = linmod('generalized_plant')
M = ltisys(Am,Bm,Cm,Dm);
M = sderiv(M,2,[1/p3 1]);
M = sderiv(M,2,[1/p3 1]);

[gopt,Gcmod1] = hinflmi(M,[1 1], 0, 0.01,[0 0 0]);

[A,B,C,D] = ltiss(Gcmod1);
Gcmod = zpk(tf(ss(A,B,C,D)))

%%

Gc = minreal(Gcmod*(1+s/2.899e06)*(s^2 + 9.848*s + 24.79)/(s*(s+0.1944)),1e-3)

%% Ln
close 6

omega_L = logspace(-5,5,10000);
L = Gcmod*Ga*Gpn*Gs*Gf;
Ln = Gc*Ga*Gpn*Gs*Gf;

figure(6)
hold on, grid on,
myngridst(Tpo,Spo)
nichols(omega_L,L,'b')
nichols(omega_L,Ln,'r')

Sn = minreal(1/(1+Ln),1e-3);
figure(1)
bodemag(omega,Sn,'r')

Tn = minreal(Ln/(1+Ln),1e-3)
figure(2)
bodemag(omega,Tn,'r')

%% Time simulation
close 7
t = linspace(0,0.015,5000);

yr = step(Tn*Kd,t);

figure(7)
hold on,grid on
plot(t,yr,'b','LineWidth',1.5);
yline(Kd*1.05);
yline(Kd*0.95);
yline((1.12)*Kd);
xline(tr,'--')
xline(ts,'-')

%%
close 8
t = linspace(0,10,10000);

da = Da0*ones(length(t),1);

yda = lsim(Gp*Sn,da,t);

figure(8)
hold on,grid on
plot(t,yda,'b','LineWidth',1.5)
yline(1.75e-2)

%%
close 9
t = linspace(0,0.1,10000);

dp = Dp0*t;

ydp = lsim(Sn,dp,t);

figure(9)
hold on,grid on
plot(t,yda,'b','LineWidth',1.5)

%%
close 10

ds = as*sin(3000*t);

yds = lsim(Tn/Gs,ds,t);

figure(10)
hold on,grid on
plot(t,yds,'b','LineWidth',1.5)
yline(2e-4)
yline(-2e-4)

%%
close 11

yu = step(Ga*Gpn*Sn,t);

figure(11)
hold on,grid on
plot(t,yu,'b','LineWidth',1.5)

%% 
figure(5)
bodemag(omega,Wu)
bodemag(omega,1/Tn)

%% mu analysis RP4
close 12
WT_mu_RP4 = tf(1/(M_HF));
omega_RP4 = logspace(log10(3000),5,1000);

[A,B,C,D] = linmod('N_scheme_RP4');
N = pck(A,B,C,D);

Nf = frsp(N,omega_RP4);

deltaset = [-1,1;-3,0;1,1]

mubounds = mu(Nf,deltaset)
figure(12)
vplot('liv,m',mubounds)

%% RP 4 U 
close 21
WT_mu_RP4 = tf(1/(M_HF));
omega_RP4 = logspace(log10(3000),5,1000);

[A,B,C,D] = linmod('N_scheme_RP4_U');
N = pck(A,B,C,D);

Nf = frsp(N,omega_RP4);

deltaset = [1,1;1,1]

mubounds = mu(Nf,deltaset)
figure(21)
vplot('liv,m',mubounds)

%% RS
K_min = 20;
K_max = 60;
K_n = (K_max + K_min)/2
Wk = (K_max - K_min)/2

b_min = 48;
b_max = 120;
b_n = (b_max + b_min)/2
Wb = (b_max - b_min)/2

c = 0.707
s = tf('s')
Gpn = zpk(K_n/(s*(1+2*c*s/b_n + (s/b_n)^2)))

omega = logspace(-2,5,2000);
Gpn_mag = squeeze(freqresp(Gpn,omega));
mag_max = zeros(length(omega),1);

figure(15)
hold on,grid on

for K = linspace(K_min,K_max,30)
    for b = linspace(b_min,b_max,30)
        Gp = K/(s*(1+2*c*s/b + (s/b)^2));
        Gp_mag = squeeze(freqresp(Gp,omega));
        delta = abs(Gp_mag./Gpn_mag -1);

        mag_max = max(delta,mag_max);
        plot(omega,delta)

    end
end

plot(omega, mag_max,'r','LineWidth',1.5)
set(gca,'XScale','log')
set(gca, 'YScale','log')

%%
figure(16)
magg = vpck(mag_max,omega);
wu = fitmag(magg);
[A,B,C,D] = unpck(wu);

Wu_new = zpk(tf(ss(A,B,C,D)))

%% RS

[A,B,C,D] = linmod('N_scheme_RS');
N = pck(A,B,C,D);

Nf = frsp(N,omega);

deltaset = [-1,1;-3,0]

mubounds = mu(Nf,deltaset)
figure(17)
vplot('liv,m',mubounds)

%% Robust stability U
[A,B,C,D] = linmod('N_scheme_RS_U');
N = pck(A,B,C,D);

Nf = frsp(N,omega);

deltaset = [1,1];

mubounds = mu(Nf,deltaset)
figure(17)
vplot('liv,m',mubounds)
%%
figure(18)
bodemag(omega,Wu_new*Tn,'--b')
hold on,grid on
bodemag(omega,Sn/WS_inv,'--r')

bodemag(omega,Wu*Tn + Sn/WS_inv,'b')

