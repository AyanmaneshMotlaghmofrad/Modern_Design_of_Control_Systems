clear variables
close all
clc
format compact

%%
z_min = 10;
z_max = 30;
z_n = (z_max+z_min)/2
Wz =  (z_max-z_min)/2

p_min = 160;
p_max = 240;
p_n = (p_max+p_min)/2
Wp = (p_max-p_min)/2

s= tf('s')
Gpn = 2e5*(s+z_n)/((s-20)*(s+p_n)*(1+s/200))

kp = dcgain(Gpn)

Gs = 5;
Ga = 1;
da0 = 8.5e-3;
dp0 = 3e-3;

as = 1e-2;

%%
Kd = 3;
Gf = 1/(Gs*Kd);

nu = 1;

M_HF = 2e-4*Gs/as
M_HF_dB = mag2db(M_HF)

%% transient requirements
s_hat = 0.08;

zeta = abs(log(s_hat))/sqrt(pi^2 + log(s_hat)^2 )
Tpo = 1/(2*zeta*sqrt(1-zeta^2))
Spo =(2*zeta*sqrt(2+4*zeta^2 + 2*sqrt(1+8*zeta^2))/(sqrt(1+8*zeta^2) + 4*zeta^2 -1))

Tpo_dB = mag2db(Tpo)

tr = 0.0625;
wn_tr = (pi-acos(zeta))/(tr*sqrt(1-zeta^2))

ts = 0.25;
wn_ts = -log(0.05)/(zeta*ts)

wn = max(wn_ts,wn_tr)
%% W1
% close 1
S2 = s*(s+2*zeta*wn)/(s^2 + 2*zeta*wn*s + wn^2)

S_star_0 = dcgain(S2/s);
omega = logspace(-3,4,2000);

figure(1)
bodemag(omega,S2,'--k')
hold on,grid on

bodemag(omega,tf(Spo),'k')

z1 = 26;
p1 = sqrt(Spo*z1/(S_star_0))
WS_inv = s*S_star_0*(1+s/z1)/(1+(2*0.82*s/p1) + (s/p1)^2)
bodemag(omega,WS_inv,'b')

z1 = 35
p1 = sqrt(Spo*z1/(S_star_0))

WS_inv1 = s*S_star_0*(1+s/z1)/(1+(2*0.707*s/p1) + (s/p1)^2)
bodemag(omega,WS_inv1,'g')

%% W2
% close 2
T2 = wn^2 /(s^2 + 2*zeta*wn*s + wn^2);

figure(2)
hold on,grid on
bodemag(omega,T2,'--k')
bodemag({2000,10000},tf(M_HF),'k')

p2 = 2000*10^((Tpo_dB-M_HF_dB)/-40);

WT_inv = Tpo/(1+ (2*0.72*s/p2) + (s/p2)^2)

bodemag(omega,WT_inv,'b')

%% Wu

Gpn_mag = squeeze(freqresp(Gpn,omega));

mag_max = zeros(length(omega),1);

figure(3)
hold on,grid on

for z = linspace(z_min,z_max,30)
    for p = linspace(p_min,p_max,30)

        Gp = 2e5*(1+s/z)/((s-20)*(s+p)*(1+s/200));

        Gp_mag = squeeze(freqresp(Gp,omega));

        delta = abs(Gp_mag./Gpn_mag -1);

        mag_max = max(delta,mag_max);

        plot(omega,delta)
    end
end

plot(omega,mag_max,'r','LineWidth',1.2)
set(gca,'XScale','log')
set(gca,'YScale','log')

%%
magg = vpck(mag_max,omega)
figure(4)
wu = fitmag(magg)
[A,B,C,D] = unpck(wu);

Wu = zpk(tf(ss(A,B,C,D)))
%%
% close 5
figure(5)
bodemag(omega,1/WT_inv)
hold on,grid on
bodemag(omega,Wu)

%% W1 and W2
W1 = 1/WS_inv1;

lambda = 0.5;
W1mod = minreal(W1*s/(s+lambda),1e-3);

figure(1)
bodemag(omega,1/W1mod,'c')

W2mod = tf(1/Tpo);
figure(2)

%%
p2 = 2300*10^((Tpo_dB-M_HF_dB)/-40);

[A,B,C,D]  = linmod('generalized_plant')
M = ltisys(A,B,C,D)
M = sderiv(M,2,[1/p2 1]);
M = sderiv(M,2,[1/p2 1]);

[gopt,Gcmod1] = hinflmi(M,[1,1],0,0.01,[0 0 0]);

[A,B,C,D] = ltiss(Gcmod1);
Gcmod = zpk(tf(ss(A,B,C,D)))

%% clearning
Gc = minreal(Gcmod*((s+0.4955))/(s),1e-3)

%%
close 6
L = Gcmod*Ga*Gpn*Gs*Gf;
Ln = Gc*Ga*Gpn*Gf*Gs;

omega_L = logspace(-5,5,10000);

figure(6)
hold on,grid on
myngridst(Tpo,Spo)
nichols(L,omega_L,'b')
nichols(Ln,omega_L)

Tn = minreal(Ln/(1+Ln),1e-3);
Sn = minreal(1/(1+Ln),1e-3);

%%
figure(1)
bodemag(omega,Sn,'r')

figure(2)
bodemag(omega,Tn,'r')

figure(8)
t = linspace(0,1,1000)
step(Tn*Kd,t)
yline(Kd*1.08)