clc
close all
clear variables
format compact
%%
K_min = 0.8;
K_max = 1.2;
K_n = (K_max+K_min)/2;
Wk =  (K_max-K_min)/2;

b_min = 420;
b_max = 580;
b_n = (b_max+b_min)/2;
Wb =  (b_max-b_min)/2;

s = tf('s');
Gpn = 10e6*K_n/(s^2 + b_n*s)
kp = dcgain(s*Gpn)

Gs = 1;
Ga = 3;

%% steady state requirements
Kd = 2
Gf = 1/(Gs*Kd)
Da0 = 8.5e-3
S_star_0 = abs(25e-3/(Da0*kp));

ap = 6e-2;
as = 1e-3;

M_LF =6e-4/ap
M_LF_dB = mag2db(M_LF)

M_HF = (Gs*1e-4)/as
M_HF_dB = mag2db(M_HF)

nu = 0;
%% transient requirements
s_hat = 0.08;
zeta = -log(s_hat)/sqrt(pi^2+log(s_hat)^2)

Spo = (2*zeta*sqrt(2+4*zeta^2 + 2*sqrt(1+8*zeta^2)))/(sqrt(1+8*zeta^2) + 4*zeta^2 -1)
Tpo = 1/(2*zeta*sqrt(1-zeta^2))
Tpo_dB = mag2db(Tpo)

tr = 5e-3
wn_tr = (pi-acos(zeta))/(tr*sqrt(1-zeta^2))

ts = 8.5e-3;
wn_ts = -log(0.05)/(ts*zeta)

wn = max([wn_tr,wn_ts])

%% W1

omega = logspace(0,5,1000);
S2 = s*(s+2*zeta*wn)/(s^2+2*zeta*wn*s+wn^2);

figure(1)
bodemag(omega,S2,'--k')
hold on, grid on
bodemag(omega,tf(Spo),'k')
bodemag({0.1,10},tf(M_LF),'k')

%WS_inv for performance
z1 = 1.45;
p1 = 23;
z2 = 350;
p2 =sqrt(Spo*z1*z2/(S_star_0*p1))
WS_inv = s*S_star_0*(1+s/z1)*(1+s/z2)/((1+s/p1)*(1+(2*0.82*s/p2)+(s/p2)^2));

bodemag(omega,WS_inv,'b')


%WS_inv for design W1
z1 = 1.45;
p1 = 20;
z2 = 350;
p2 =sqrt(Spo*z1*z2/(S_star_0*p1))
WS_inv1 = s*S_star_0*(1+s/z1)*(1+s/z2)/((1+s/p1)*(1+(2*0.707*s/p2)+(s/p2)^2));

bodemag(omega,WS_inv1,'g')

%% W2

T2 = wn^2/(s^2 + 2*zeta*wn*s + wn^2);

figure(2)
bodemag(omega,T2,'--k')
hold on, grid on
bodemag({10000,100000},tf(M_HF),'k')
bodemag(omega, tf(Tpo))

% WT_inv
p3 = 10000*10^((Tpo_dB-M_HF_dB)/-40);
WT_inv = Tpo/(1+(2*s*0.82/p3)+(s/p3)^2)

bodemag(omega,WT_inv,'b')


%% Wu
 
Gpn_mag = squeeze(freqresp(Gpn,omega));
mag_max = zeros(length(omega),1);

figure(3)
hold on, grid on,

for K = linspace(K_min,K_max,30)
    for b = linspace(b_min,b_max,30)
        Gp = 10e6*K/(s^2+b*s);
        Gp_mag = squeeze(freqresp(Gp,omega));

        delta = abs(Gp_mag./Gpn_mag-1);
        plot(omega,delta)
        mag_max = max(delta,mag_max);
    end
end

plot(omega,mag_max,'r','LineWidth',1.5)
set(gca,'XScale','log')
set(gca,'YScale','log')

%% Wfigure,
magg = vpck(mag_max,omega)
W_u_m = fitmag(magg);
[A,B,C,D] = unpck(W_u_m);
Wu = zpk(tf(ss(A,B,C,D)));

%%

figure(5)
bodemag(omega, Wu,'r')
hold on, grid on
bodemag(omega, WT_inv,'b')

%% W1 W2

W2 = (1+s/p3)^2/Tpo
figure(2)
bodemag(omega, 1/W2,'g')

W2mod = tf(1/Tpo);


W1 = 1/WS_inv1;
lambda = 4;
W1mod = minreal(W1*s/(s+lambda),1e-4);

figure(1)
bodemag(omega,1/W1mod,'c')
%% Hinf
[Am,Bm,Cm,Dm] = linmod('generalized_plant')
M = ltisys(Am,Bm,Cm,Dm);
M = sderiv(M,2,[1/p3 1]);
M = sderiv(M,2,[1/p3 1]);
[gopt,Gcmod1] = hinflmi(M,[1 1],0,0.01,[0 0 0]);

[A,B,C,D] = ltiss(Gcmod1);
Gcmod = zpk(tf(ss(A,B,C,D)))
%% clearning the controller
Gc = minreal(Gcmod*(1+s/8.482e06)*(s+2.933)/((s+0.1102)),1e-3)


%%
L = Gcmod *Ga*Gpn*Gs*Gf;
omega_L= logspace(-2,5,10000);
figure(6)
nichols(L,omega)
hold on,grid on

Ln = Gc*Ga*Gpn*Gs*Gf
nichols(Ln,omega)

myngridst(Tpo,Spo)

Tn = minreal(Ln/(1+Ln),1e-3)
figure(2)
bodemag(omega,Tn,'r')

Sn = minreal(1/(1+Ln),1e-3)
figure(1)
bodemag(omega,Sn,'r')
%% 
close 7
t = linspace(0,0.2,1000);
yr = step(Tn*Kd,t);

figure(7)
plot(t,yr,'b','LineWidth',1.5)
hold on,grid on
yline(Kd*1.05)
yline(Kd*0.95)
xline(5e-3,'b')
xline(8.5e-3,'r')

%% error
da = Da0*ones(length(t),1);

yda = lsim(Gpn*Sn,da,t);

figure(8)
hold on, grid on
plot(t,yda,'b','LineWidth',1.5)
yline(25e-3)


%% dp
tdp = linspace(0,10,2000);
dp = ap*sin(tdp*10)

ydp = lsim(Sn,dp,tdp)

figure(9)
hold on, grid on
plot(tdp,ydp,'b','LineWidth',1.5)
yline(6e-4)
yline(-6e-4)

%% ds
tds = linspace(0,0.5,1000);
ds = as*sin(tds*10000)

yds = lsim(Tn/Gs,ds,tds)

figure(10)
hold on, grid on
plot(tds,yds,'b','LineWidth',1.5)
yline(8.5e-4)
yline(-8.5e-4)

%% mu analysis
omega_RP4 = logspace(4,6,1000);
WT_mu_RP4 = tf(1/M_HF)
[An,Bn,Cn,Dn] = linmod('N_scheme_RP4')
N = pck(An,Bn,Cn,Dn);
Nf = frsp(N,omega_RP4);
deltaset=[-1,1;-1,1;1,1];

mubounds = mu(Nf,deltaset)

figure(11)
vplot('liv,m',mubounds)

%% Robust stability
K_min = 0.4;
K_max = 1.6;
K_n = (K_max+K_min)/2;
Wk =  (K_max-K_min)/2;

b_min = 400;
b_max = 600;
b_n = (b_max+b_min)/2;
Wb =  (b_max-b_min)/2;

omega_RS = logspace(-2,6,2000);
[An,Bn,Cn,Dn] = linmod('N_scheme_RS')
N = pck(An,Bn,Cn,Dn);
Nf = frsp(N,omega_RS);
deltaset=[-1,1;-1,1];

mubounds = mu(Nf,deltaset)

figure(12)
vplot('liv,m',mubounds)

%% Wunew
Gpn_mag = squeeze(freqresp(Gpn,omega));
mag_max = zeros(length(omega),1);

figure(14)
hold on, grid on,

for K = linspace(K_min,K_max,30)
    for b = linspace(b_min,b_max,30)
        Gp = 10e6*K/(s^2+b*s);
        Gp_mag = squeeze(freqresp(Gp,omega));

        delta = abs(Gp_mag./Gpn_mag-1);
        plot(omega,delta)
        mag_max = max(delta,mag_max);
    end
end

plot(omega,mag_max,'r','LineWidth',1.5)
set(gca,'XScale','log')
set(gca,'YScale','log')

%% Wfigure,
magg = vpck(mag_max,omega)
figure(15)
W_u_m = fitmag(magg);
[A,B,C,D] = unpck(W_u_m);
Wu_new = zpk(tf(ss(A,B,C,D)))

%%

figure(16)
bodemag(omega, 1/Wu_new,'r')
hold on, grid on
bodemag(omega, Tn,'b')

figure(16)
bodemag(omega,Wu*Tn,'b')
hold on, grid on,
bodemag(omega,Sn/WS_inv,'r')
bodemag(omega,Wu*Tn + Sn/WS_inv,'g')














