clc
clear variables
close all 
format compact
for i = 1:20
    figure(20)
end
%% data

K_min = 1.4;
K_max = 2.6;
K_n = (K_max+K_min)/2
Wk = (K_max-K_min)/2

p_min = 75;
p_max = 175;
p_n = (p_max+p_min)/2
Wp = (p_max-p_min)/2

s= tf('s');
Gpn = K_n/(s*(1+s/p_n))
kp = dcgain(s*Gpn)

Gs = 2
Ga = -12500
Da0 = 8.5e-2
ap = 2e-2
as = 5e-2

%% steady state requirements
Kd = 2;

Gf = 1/(Kd*Gs)

M_LF = 2e-3/ap

M_HF = 2.5e-3*Gs/as
M_HF_dB = mag2db(M_HF)

%%
s_hat = 0.1;
zeta = -log(s_hat)/sqrt(pi^2 + log(s_hat)^2)

Tpo = 1/(2*zeta*sqrt(1-zeta^2))
Tpo_dB = mag2db(Tpo)

Spo = (2*zeta*sqrt(2+4*zeta^2 + 2*sqrt(1+8*zeta^2)))/(sqrt(1+8*zeta^2)+4*zeta^2 -1)
Spo_dB = mag2db(Spo)

tr = 0.06485

wn_tr = (pi-acos(zeta))/(tr*sqrt(1-zeta^2))

ts = 0.1238

wn_ts = -log(0.05)/(zeta*ts)

wn = max(wn_tr,wn_ts)

%% WS_inv
close 1
S2 = s*(s+ 2*zeta*wn)/(s^2 + 2*zeta*wn*s + wn^2);

omega = logspace(-3,4,2000);

figure(1)
hold on,grid on
bodemag({0.01,2.5},tf(M_LF),'k')
bodemag(omega, tf(Spo),'k')
bodemag(omega,S2,'--k')

% WS for performance
S_star_0 = 1;
p1 = 0.028;
z1 = 20;
p2 = sqrt(Spo*z1/p1*S_star_0);
WS_inv = S_star_0*s^2*(1+s/z1)/((1+s/p1)*(1+2*0.82*s/p2 + (s/p2)^2));
bodemag(omega,WS_inv,'b')

% WS for design
S_star_0 = 1;
p1 = 0.027;
z1 = 30;
p2 = sqrt(Spo*z1/p1*S_star_0);
WS_inv1 = S_star_0*s^2*(1+s/z1)/((1+s/p1)*(1+2*0.72*s/p2 + (s/p2)^2));
bodemag(omega,WS_inv1,'g')

%% WT
close 2
T2 = wn^2/(s^2 + 2*zeta*wn*s + wn^2);

figure(2)
hold on,grid on
bodemag({250,1000},tf(M_HF),'k')
bodemag(omega,tf(Tpo),'k')
bodemag(omega,T2,'--k')

p3 = 250*10^((Tpo_dB - M_HF_dB)/-40)

WT_inv = Tpo/(1+ 2*0.78*s/p3 + (s/p3)^2)

bodemag(omega,WT_inv,'b')

%% 
close 3
Gpn_mag = squeeze(freqresp(Gpn,omega));

mag_max = zeros(length(omega),1);

figure(3)
hold on,grid on

for K = linspace(K_min,K_max,30)
    for p = linspace(p_min,p_max,30)
        
        Gp = K/(s*(1+s/p));

        Gp_mag = squeeze(freqresp(Gp,omega));

        delta = abs(Gp_mag./Gpn_mag -1);

        mag_max = max(delta,mag_max);

        plot(omega,delta)
    end
end
plot(omega,mag_max,'r','LineWidth',1.5)
set(gca,'XScale','log')
set(gca,'YScale','log')

%% WU
magg = vpck(mag_max,omega);

wu = fitmag(magg);
[A,B,C,D] = unpck(wu);
Wu = zpk(tf(ss(A,B,C,D)))

%% Wu and WT

figure(5)
hold on,grid on
bodemag(omega, Wu,'--b')
bodemag(omega,1/WT_inv,'b')

%% W1 and W2
W2 = (1+s/p3)^2/Tpo
W2mod = tf(1/Tpo);

figure(2)
bodemag(omega,1/W2,'c')

W1 = 1/WS_inv1;
lambda = 0.1;
W1mod = minreal(W1*s^2/(s+lambda)^2,1e-3)

figure(1)
bodemag(omega,1/W1mod,'c')

%%  Hinf
[A,B,C,D] = linmod('generalized_plant')
M = ltisys(A,B,C,D);

M = sderiv(M,2,[1/p3 1]);
M = sderiv(M,2,[1/p3 1]);

[gopt,Gcmod1] = hinflmi(M,[1 1],0,0.01,[0 0 0])
[A,B,C,D] = ltiss(Gcmod1);
Gcmod = zpk(tf(ss(A,B,C,D)))

%% clearning
Gc = minreal(Gcmod*(s^2 + 0.1793*s + 0.00949)*(1+s/2.323e05)/(s*(s+0.02707)),1e-3)

%% L
close 6
omega_L = logspace(-5,5,5000);
L = Gcmod*Ga*Gpn*Gs*Gf;
Ln = Gc*Ga*Gpn*Gs*Gf;

figure(6)
hold on,grid on
myngridst(Tpo,Spo)
nichols(L,omega_L,'b')
nichols(Ln,omega_L,'r')

Sn =minreal(1/(1+Ln),1e-3)
figure(1)
bodemag(omega,Sn,'r')

Tn = minreal(Ln/(1+Ln),1e-3)
figure(2)
bodemag(omega,Tn,'r')
 %% time simulation


t = linspace(0,0.2,1000);

yr = step(Kd*Tn,t)

figure(7)
hold on,grid on
plot(t,yr,'b','LineWidth',1.5)
yline(Kd*1.05,'-')
yline(Kd*0.95,'-')
yline(Kd,'--')
xline(tr,'--')
xline(ts,'-')
%% u
u = step(Gc*Ga*Sn,t);

figure(9)
hold on,grid on
plot(t,u,'b','LineWidth',1.5)

%% dp
close 11

ds = as*sin(250*t);

yds = lsim(Tn/Gs,ds,t);

figure(11)
hold on,grid on
plot(t,yds,'b','LineWidth',1.5)
yline(2.5e-3)
yline(-2.5e-3)

%%
t = linspace(0,10,10000)
dp = ap*sin(2.5*t);

ydp = lsim(Sn,dp,t);

figure(10)
hold on ,grid on
plot(t,ydp,'b','LineWidth',1.5)
yline(2e-3)
yline(-2e-3)
%% da 
t = linspace(0,100,10000)
da = Da0*ones(length(t),1);

yda = lsim(Gpn*Sn,da,t);

figure(8)
hold on,grid on
plot(t,yda,'b','LineWidth',1.5)

%% Robust performance
figure(12)
hold on,grid on
bodemag(omega,Wu*Tn,'--b')
bodemag(omega,Sn/WS_inv,'--r')
bodemag(omega, Wu*Tn + Sn/WS_inv,'r')

%% mu analysis

WT_mu_RP4 = tf(1/M_HF);

omega_RP4 = logspace(log10(250),5,1000);

[A,B,C,D] = linmod('N_scheme_RP4');
N = pck(A,B,C,D);

Nf = frsp(N,omega_RP4);
deltaset =[-1,1;-1,1;1,1];

mubounds = mu(Nf,deltaset);

figure(13)
vplot('liv,m',mubounds)

%%

[A,B,C,D] = linmod('N_scheme_RP4_U');
N = pck(A,B,C,D);

Nf = frsp(N,omega_RP4);
deltaset =[1,1;1,1];

mubounds = mu(Nf,deltaset);

figure(14)
vplot('liv,m',mubounds)

%% robust stability

K_min = 1;
K_max = 3;
K_n = (K_max+K_min)/2
Wk = (K_max-K_min)/2

p_min = 50;
p_max = 200;
p_n = (p_max+p_min)/2
Wp = (p_max-p_min)/2

s= tf('s');
Gpn = K_n/(s*(1+s/p_n))
kp = dcgain(s*Gpn)
[A,B,C,D] = linmod('N_scheme_RS');
N = pck(A,B,C,D);

Nf = frsp(N,omega);
deltaset =[-1,1;-1,1];

mubounds = mu(Nf,deltaset);

figure(15)
vplot('liv,m',mubounds)

%% RS U
Gpn_mag = squeeze(freqresp(Gpn,omega));

mag_max = zeros(length(omega),1);

figure(16)
hold on,grid on

for K = linspace(K_min,K_max,30)
    for p = linspace(p_min,p_max,30)
        
        Gp = K/(s*(1+s/p));

        Gp_mag = squeeze(freqresp(Gp,omega));

        delta = abs(Gp_mag./Gpn_mag -1);

        mag_max = max(delta,mag_max);

        plot(omega,delta)
    end
end
plot(omega,mag_max,'r','LineWidth',1.5)
set(gca,'XScale','log')
set(gca,'YScale','log')

%% WU
magg = vpck(mag_max,omega);

wu = fitmag(magg);
[A,B,C,D] = unpck(wu);
Wu_new = zpk(tf(ss(A,B,C,D)))

%%
[A,B,C,D] = linmod('N_scheme_RS_U');
N = pck(A,B,C,D);

Nf = frsp(N,omega);
deltaset =[1,1];

mubounds = mu(Nf,deltaset);

figure(16)
vplot('liv,m',mubounds)



%%
figure(17)
Wunew_mag = abs(squeeze(freqresp(Wu_new*Tn,omega)));

plot(omega,Wunew_mag,'b','LineWidth',1.5)
set(gca,'XScale','log')



