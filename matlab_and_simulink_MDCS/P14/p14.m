clc
clear variables
close all
format compact
for i = 1:14
    figure(i)
end

%% Given data

a_min = 110;
a_max = 142;
a_n = (a_max + a_min)/2
Wa =  (a_max - a_min)/2

b_min = 36;
b_max = 44;
b_n = (b_max + b_min)/2
Wb =  (b_max - b_min)/2

c_min = 0.095;
c_max = 0.105;
c_n = (c_max + c_min)/2
Wc =  (c_max - c_min)/2
s = tf('s');
Gpn = (b_n-1)*a_n^2/(a_n^2+ 14*a_n*c_n*s + s^2)


kp = dcgain(Gpn)

Gs = 1;
Ga = 0.01;
Da0 = 8.5e-3;
Dp0 = 3e-3;
as = 1e-2;

%% steady state requirements
Kd = 6;
Gf = 1/(Kd*Gs)

S_star_0 = (1.75e-2)/(kp*Da0) %nu 0



nu = 1;
M_HF = (1e-4*Gs)/(as)
M_HF_dB = mag2db(M_HF)

%% transient requirements
s_hat = 0.09;

zeta = -log(s_hat)/sqrt(pi^2+log(s_hat)^2)
Spo = (2*zeta*sqrt(2+4*zeta^2+2*sqrt(1+8*zeta^2)))/(sqrt(1+8*zeta^2) + 4*zeta^2 - 1)
Tpo = 1/(2*zeta*sqrt(1-zeta^2))
Tpo_dB = mag2db(Tpo);

tr =0.0065;
wn_tr = (pi-acos(zeta))/(tr*sqrt(1-zeta^2))

ts =0.0105;
wn_ts = -log(0.05)/(ts*zeta)

wn = max([wn_tr,wn_ts])

%% W1
close 1
S2 = s*(s+2*zeta*wn)/(s^2+2*zeta*wn*s+wn^2);

S_star_0 = dcgain(S2/s);
omega = logspace(0,5,1000);

figure(1)
hold on , grid on
bodemag(omega,S2,'--k')
bodemag(omega,tf(Spo),'k')

% WS_inv for the performance
z1 = 300;
p1 = sqrt(Spo*z1/S_star_0)
WS_inv = S_star_0*s*(1+s/z1)/(1+ (2*0.8*s/p1)+ (s/p1)^2)
bodemag(omega,WS_inv,'b')

% WS_inv1 for W1
z1 = 350;
p1 = sqrt(Spo*z1/S_star_0)
WS_inv1 = S_star_0*s*(1+s/z1)/(1+ (2*0.707*s/p1)+ (s/p1)^2)
bodemag(omega,WS_inv1,'g')

%% WT_inv
close 2
T2 = wn^2/(s^2+2*zeta*wn*s+wn^2);
p2 = 12000*10^((Tpo_dB-M_HF_dB)/-40);

figure(2)
bodemag({12000,100000},tf(M_HF),'k')
hold on, grid on
bodemag(omega,T2,'--k')
bodemag(omega,tf(Tpo),'k')


WT_inv = zpk(Tpo/(1+(2*0.707*s/p2)+ (s/p2)^2))
bodemag(omega,WT_inv,'b')

%% Wu

Gpn_mag = squeeze(freqresp(Gpn,omega));
mag_max = zeros(length(omega),1);

figure(3)
hold on, grid on

for b = linspace(b_min,b_max,10)
    for a = linspace(a_min,a_max,10)
        for c = linspace(c_min,c_max,10)
            Gp = (b-1)*a^2/(a^2+14*a*c*s+s^2);

            Gp_mag = squeeze(freqresp(Gp,omega));
            
            delta = abs(Gp_mag./Gpn_mag - 1);
            mag_max = max(delta,mag_max);
            plot(omega,delta)
        end
    end
end

plot(omega,mag_max,'r','LineWidth',1.5)
set(gca,'XScale','log')
set(gca,'YScale','log')

%%
magg = vpck(mag_max,omega);
figure(4)
wu = fitmag(magg);
[A,B,C,D] = unpck(wu);

Wu = zpk(tf(ss(A,B,C,D)))

%% Wu and WT_inv
close 5
figure(5)
bodemag(omega,Wu,'b')
hold on, grid on
bodemag(omega,1/WT_inv,'r')


%% W1 and W2

W2 = (1+s/p2)^2/Tpo;
W2mod = tf(1/Tpo)

% W1
W1 = 1/WS_inv1;
lambda = 1;

W1mod = minreal(zpk(W1*s/(s+lambda)),1e-4)
figure(1)
bodemag(omega,1/W1mod,'c')

%% Hinf
[Am,Bm,Cm,Dm] = linmod('generalized_plant')
M = ltisys(Am,Bm,Cm,Dm);
M = sderiv(M,2,[1/p2 1]);
M = sderiv(M,2,[1/p2 1]);

[gopt,Gcmod1] = hinflmi(M,[1 1],0,0.01,[0 0 0]);
[A,B,C,D] = ltiss(Gcmod1);

Gcmod = zpk(tf(ss(A,B,C,D)))

%% cleaning
Gc = minreal(Gcmod*(s+0.6977)*(1+s/7.901e06)/s,1e-4)
%%
close 6
omega_L = logspace(-4,6,50000);
L = Gcmod*Ga*Gpn*Gs*Gf;
Ln = Gc*Ga*Gpn*Gs*Gf;

figure(6)
hold on,grid on
myngridst(Tpo,Spo)
nichols(L,omega_L)
nichols(Ln,omega_L)

%% Sn
Sn = minreal(1/(1+Ln),1e-3);

figure(1)
bodemag(omega,Sn,'r')
%% Tn
Tn = minreal(Ln/(1+Ln),1e-4);

figure(2)
bodemag(omega,Tn,'r')

%% time domain simulation
t = linspace(0,0.15,2000);
yr = step(Kd*Tn,t);



figure(7)
hold on ,grid on
plot(t,yr,'b','LineWidth',1.5)
xline(tr,'--r')
xline(ts,'r')
yline(Kd*1.05)
yline(Kd*0.95)
yline(1.09*Kd)

%% da
da = Da0*ones(length(t),1);
yda = lsim(Gp*Sn,da,t);

figure(8)
hold on,grid on
plot(t,yda,'b','LineWidth',1.5)
yline(1.75e-2)

%%dp

dp = Dp0*t;
ydp = lsim(Sn,dp,t);

figure(9)
hold on,grid on,
plot(t,ydp,'b','LineWidth',1.5)

%% ds
ds =as*sin(12000*t);
yds = lsim(Tn/Gs,ds,t);
figure(10)
hold on, grid on
plot(t,yda,'b','LineWidth',1.5)
yline(1e-4)
yline(-1e-4)

%% Robust performance

WT_mu_RP4 = tf(1/M_HF);
omega_RP4 = logspace(log10(12000),6,1000);

[An,Bn,Cn,Dn] = linmod('N_scheme_RP4');
N = pck(An,Bn,Cn,Dn);
Nf = frsp(N,omega_RP4);
deltaset = [-1,1;-5,0;-1,1;1,1];

mubounds = mu(Nf,deltaset);

vplot('liv,m',mubounds)

%% RS

a_min = 104;
a_max = 148;
a_n = (a_max + a_min)/2
Wa =  (a_max - a_min)/2

b_min = 28;
b_max = 52;
b_n = (b_max + b_min)/2
Wb =  (b_max - b_min)/2

c_min = 0.085;
c_max = 0.115;
c_n = (c_max + c_min)/2
Wc =  (c_max - c_min)/2
omega_RS = logspace(-2,6,2000);

[An,Bn,Cn,Dn] = linmod('N_scheme_RS');
N = pck(An,Bn,Cn,Dn);
Nf = frsp(N,omega_RS);
deltaset = [-1,1;-5,0;-1,1];

mubounds = mu(Nf,deltaset);

figure,
vplot('liv,m',mubounds)

%% Wunew

Gpn_mag = squeeze(freqresp(Gpn,omega));
mag_max = zeros(length(omega),1);

figure(3)
hold on, grid on

for b = linspace(b_min,b_max,10)
    for a = linspace(a_min,a_max,10)
        for c = linspace(c_min,c_max,10)
            Gp = (b-1)*a^2/(a^2+14*a*c*s+s^2);

            Gp_mag = squeeze(freqresp(Gp,omega));
            
            delta = abs(Gp_mag./Gpn_mag - 1);
            mag_max = max(delta,mag_max);
            plot(omega,delta)
        end
    end
end

plot(omega,mag_max,'r','LineWidth',1.5)
set(gca,'XScale','log')
set(gca,'YScale','log')

%% Wu_new
magg = vpck(mag_max,omega);
figure(4)
wu = fitmag(magg);
[A,B,C,D] = unpck(wu);

Wu_new = zpk(tf(ss(A,B,C,D)))


%% RP
figure(13)
bodemag(omega,Wu_new*Tn,'g')
grid on, hold on
bodemag(omega,Sn/WS_inv,'r')
bodemag(omega,Wu_new*Tn + Sn/WS_inv,'b')


omega_RS = logspace(-2,6,2000);

[An,Bn,Cn,Dn] = linmod('N_scheme_RS_unstructured');
N = pck(An,Bn,Cn,Dn);
Nf = frsp(N,omega_RS);
deltaset = [1,1];

mubounds = mu(Nf,deltaset);

figure,
vplot('liv,m',mubounds)