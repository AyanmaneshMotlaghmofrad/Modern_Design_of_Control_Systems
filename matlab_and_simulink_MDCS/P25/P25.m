clc
clear variables
close all
format compact

for i = 1:20
    figure(i)
end

%% Data

K_min = 24;
K_max = 56;
K_n = (K_max + K_min)/2
Wk = (K_max - K_min)/2

b_min = 68;
b_max = 100;
b_n = (b_max + b_min)/2
Wb = (b_max - b_min)/2
s = tf('s');
Gpn = K_n/(1+1.414*s/b_n + (s/b_n)^2)
kp = dcgain(Gpn)

Gs = 1;
Ga = 100;
Da0 = 8.5e-3;
Dpo = 3.5e-3;
as = 1e-2;

%% steady state requirements
Kd = 1;
Gf = 1/(Kd*Gs);

S_star_0 = 3.5e-6/3.5e-3

M_HF = 2e-4*Gs/(as)
M_HF_dB = mag2db(M_HF)

%% transient
s_hat = 0.09;

zeta = -log(s_hat)/sqrt(pi^2 + log(s_hat)^2)
Tpo = 1/(2*zeta*sqrt(1-zeta^2))
Tpo_dB = mag2db(Tpo)
Spo = (2*zeta*sqrt(2+4*zeta^2 + 2*sqrt(1+8*zeta^2)))/(sqrt(1+8*zeta^2)+ 4*zeta^2 -1)

tr = 0.012;

wn_tr = (pi-acos(zeta))/(tr*sqrt(1-zeta^2))

ts = 0.0211;

wn_ts = -log(0.05)/(ts*zeta)

wn = max(wn_ts,wn_tr)

%% W1
close 1
S2 = s*(s+ 2*zeta*wn)/(s^2 + 2*zeta*wn*s + wn^2);
omega = logspace(-2,5,2000);

figure(1)
hold on, grid on
bodemag(omega,S2,'--k')
bodemag(omega,tf(Spo),'k')

%WS_inv for performance 

z1 = 0.19
z2 = 150
p1 = sqrt(Spo*z1*z2/S_star_0)
WS_inv = S_star_0*(1+s/z1)*(1+s/z2)/(1+(2*0.75*s/p1)+ (s/p1)^2);

bodemag(omega,WS_inv,'b')

% WS_inv for W1
z1 = 0.195
z2 = 180
p1 = sqrt(Spo*z1*z2/S_star_0)
WS_inv1 = S_star_0*(1+s/z1)*(1+s/z2)/(1+(2*0.707*s/p1)+ (s/p1)^2);

bodemag(omega,WS_inv1,'g')

%%
close 2

T2 = wn^2/(s^2 + 2*zeta*wn*s + wn^2);

figure(2)
hold on, grid on
bodemag(omega,T2,'--k')
bodemag({3000,10000},tf(M_HF),'k')
bodemag(omega,tf(Tpo),'k')

p2 = 3000*10^((Tpo_dB - M_HF_dB)/-40)
WT_inv = Tpo/(1 + (2*0.707*s/p2)+(s/p2)^2);

bodemag(omega,WT_inv,'b')

%% Wu
close 3
Gpn_mag = squeeze(freqresp(Gpn,omega));

mag_max = zeros(length(omega),1);

figure(3)
hold on, grid on

for K = linspace(K_min,K_max,30)
    for b = linspace(b_min,b_max,30)
        
        Gp = K/(1+(1.414*s/b) + (s/b)^2);
        Gp_mag = squeeze(freqresp(Gp,omega));

        delta = abs(Gp_mag./Gpn_mag -1);

        mag_max = max(delta, mag_max);

        plot(omega,delta)

    end
end
plot(omega,mag_max,'b','LineWidth',1.5)
set(gca, 'XScale','log')
set(gca,'YScale','log')

%%
magg = vpck(mag_max,omega);
figure(4)
wu = fitmag(magg);
[A,B,C,D] = unpck(wu);

Wu = zpk(tf(ss(A,B,C,D)))

%% Wu and WT

figure(5)
hold on,grid on
bodemag(omega, Wu,'b')
bodemag(omega, 1/WT_inv,'r')

%% W1 and W2
W2 = (1+s/p2)^2/Tpo;
figure(2)
bodemag(omega,1/W2,'c')

W2mod = tf(1/Tpo);
figure(1)
W1 = 1/WS_inv1;
bodemag(omega,1/W1,'c')

%% design
[A,B,C,D] = linmod('generalized_plant');
M = ltisys(A,B,C,D);
M = sderiv(M,2,[1/p2 1]);
M = sderiv(M,2,[1/p2 1]);

[gopt,Gcmod1] = hinflmi(M,[1 1], 0 , 0.01, [0 0 0]);

[A,B,C,D] = ltiss(Gcmod1);

Gcmod = zpk(tf(ss(A,B,C,D)))

%% clearning

Gc = minreal(Gcmod*(1+s/1.043e06),1e-3)

%%
L = Gcmod*Ga*Gpn*Gs*Gf;

Ln = Gc*Ga*Gpn*Gs*Gf;

omega_L = logspace(-5,5,10000);

figure(6)
hold on ,grid on
myngridst(Tpo,Spo)
nichols(L,omega_L,'b')
nichols(Ln,omega_L,'r')

Tn = minreal(Ln/(1+Ln),1e-3);
figure(2)
bodemag(omega,Tn,'r')

Sn = minreal(1/(1+Ln),1e-3);
figure(1)
bodemag(omega,Sn,'r')

figure(5)
bodemag(omega,1/Tn,'c')
%%
close 7
figure(7)
hold on,grid on
bodemag(omega,Wu*Tn + Sn/WS_inv,'b')
bodemag(omega,Wu*Tn,'--b')
bodemag(omega,Sn/WS_inv,'--r')
%% time simulation
%% yr
close 8

t = linspace(0,0.1,5000);

yr = step(Kd*Tn,t);

figure(8)
hold on,grid on
plot(t,yr,'b','LineWidth',1.5)
yline(1.05)
yline(0.95)
xline(0.012,'--')
xline(0.0211,'-')

%% dp
close 9

dp = Dpo*ones(length(t),1);
ydp = lsim(Sn,dp,t);

figure(9)
hold on ,grid on

plot(t,ydp,'b','LineWidth',1.5);
yline(3.5e-6)
yline(-3.5e-6)

%% ds
close 10

ds = as*sin(3000*t);

yds = lsim(Tn/Gs,ds,t);

figure(10)
hold on,grid on
plot(t,yds,'b','LineWidth',1.5)

yline(2e-4)
yline(-2e-4)

%% u
close 11
yu = step(Gpn*Ga*Sn,t);

figure(11)
hold on,grid on
plot(t,yu,'b','LineWidth',1.5)

%% mu analysis
close 12
WT_mu_RP4 = tf(1/M_HF);
omega_RP4 = logspace(log10(3000),5,1000)

[A,B,C,D] = linmod('N_scheme_RP4')
N = pck(A,B,C,D);
Nf = frsp(N,omega_RP4);
deltaset = [-1,1;-3,0;1,1];

mubounds = mu(Nf,deltaset)

figure(12)
vplot('liv,m',mubounds);

%% RP4 U
close 15
WT_mu_RP4 = tf(1/M_HF);
omega_RP4 = logspace(log10(3000),5,1000)

[A,B,C,D] = linmod('N_scheme_RP4_U')
N = pck(A,B,C,D);
Nf = frsp(N,omega_RP4);
deltaset = [1,1;1,1];

mubounds = mu(Nf,deltaset)

figure(15)
vplot('liv,m',mubounds);

%% RS
close 13
K_min = 20;
K_max = 60;
K_n = (K_max + K_min)/2
Wk = (K_max - K_min)/2

b_min = 110;
b_max = 160;
b_n = (b_max + b_min)/2
Wb = (b_max - b_min)/2
s = tf('s');
Gpn = K_n/(1+1.414*s/b_n + (s/b_n)^2)
kp = dcgain(Gpn)


[A,B,C,D] = linmod('N_scheme_RS')
N = pck(A,B,C,D);
Nf = frsp(N,omega);
deltaset = [-1,1;-3,0];

mubounds = mu(Nf,deltaset)

figure(13)
vplot('liv,m',mubounds);

%% 

close 16
Gpn_mag = squeeze(freqresp(Gpn,omega));

mag_max = zeros(length(omega),1);

figure(16)
hold on, grid on

for K = linspace(K_min,K_max,30)
    for b = linspace(b_min,b_max,30)
        
        Gp = K/(1+(1.414*s/b) + (s/b)^2);
        Gp_mag = squeeze(freqresp(Gp,omega));

        delta = abs(Gp_mag./Gpn_mag -1);

        mag_max = max(delta, mag_max);

        plot(omega,delta)

    end
end
plot(omega,mag_max,'b','LineWidth',1.5)
set(gca, 'XScale','log')
set(gca,'YScale','log')

%%
magg = vpck(mag_max,omega);
figure(17)
wu = fitmag(magg);
[A,B,C,D] = unpck(wu);

Wu_new = zpk(tf(ss(A,B,C,D)))

%% RS U with Wunew
close 17
[A,B,C,D] = linmod('N_scheme_RS_U')
N = pck(A,B,C,D);
Nf = frsp(N,omega);
deltaset = [1,1];

mubounds = mu(Nf,deltaset)

figure(17)
vplot('liv,m',mubounds);


%%
close 18
figure(18)
hold on
grid on
magnitude = bodemag(omega,Wu_new*Tn,'b')







