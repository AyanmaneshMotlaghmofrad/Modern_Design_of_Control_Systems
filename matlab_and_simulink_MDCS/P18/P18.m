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
WK = (K_max-K_min)/2

b_min = 68;
b_max = 100;
b_n = (b_max + b_min)/2
Wb = (b_max-b_min)/2
c = 0.707;
s = tf('s')
Gpn = K_n/(s*(1+2*c*s/b_n +(s/b_n)^2))
kp = dcgain(s*Gpn)
Gs = 2
Ga = 0.01;
Da0 = 8.5e-3
Dp0 = 3.5e-3
as = 1e-2

%% steady state requirements
Kd = 1.5;
Gf = 1/(Kd*Gs)

S_star_0 = 1.75e-2/(Da0*kp)

M_HF = 2e-4*Gs/as
M_HF_dB = mag2db(M_HF)

%% transient requirements

s_hat = 0.08;
zeta = -log(s_hat)/sqrt(pi^2 + log(s_hat)^2)

Tpo = 1/(2*zeta*sqrt(1-zeta^2))
Tpo_dB = mag2db(Tpo)

Spo = (2*zeta*sqrt(2+4*zeta^2 + 2*sqrt(1+8*zeta^2)))/(sqrt(1+8*zeta^2) + 4*zeta^2 -1)


tr = 0.01;
wn_tr = (pi - acos(zeta))/(tr*sqrt(1-zeta^2))

ts = 0.017;
wn_ts = -log(0.05)/(ts*zeta)

wn = max(wn_tr,wn_ts)

%% WS_inv 
% close 1

S2 = s*(s+2*zeta*wn)/(s^2 + 2*zeta*wn*s + wn^2)


omega = logspace(-3,4,2000);

figure(1)
hold on,grid on
bodemag(omega,tf(Spo),'k')
bodemag(omega,S2,'--k')

%WS for performance
p1 = 0.085;
z1 = 150;
p2 = sqrt(Spo*z1/(p1*S_star_0))
WS_inv = S_star_0*s^2*(1+s/z1)/((1+s/p1)*(1+ (2*0.8*s/p2)+(s/p2)^2))

bodemag(omega,WS_inv,'b')

% Ws for desing
p1 = 0.085;
z1 = 310;
p2 = sqrt(Spo*z1/(p1*S_star_0*0.75))
WS_inv1 = S_star_0*0.75*s^2*(1+s/z1)/((1+s/p1)*(1+ (2*0.707*s/p2)+(s/p2)^2))

bodemag(omega,WS_inv1,'g')

%% WT
close 2
T2 = wn^2 /(s^2 + 2*zeta*wn*s + wn^2)
p3 = 2000*10^((Tpo_dB-M_HF_dB)/-40)
WT_inv = Tpo/(1+ 2*0.71*s/p3 + (s/p3)^2)

figure(2)
hold on,grid on
bodemag({2000,10000},tf(M_HF),'k')
bodemag(omega,tf(Tpo),'k')
bodemag(omega,T2,'--k')
bodemag(WT_inv)

%% Wu
close 3
Gpn_mag = squeeze(freqresp(Gpn,omega));

mag_max = zeros(length(omega),1);

figure(3)
hold on,grid on
for K = linspace(K_min,K_max,30)
    for b =linspace(b_min,b_max,30)
        
        Gp = K/(s*(1+ 2*c*s/b + (s/b)^2));
        Gp_mag = squeeze(freqresp(Gp,omega));
        
        delta = abs(Gp_mag./Gpn_mag -1);

        mag_max = max(delta,mag_max);

        plot(omega,delta)
    
    
    end
end

plot(omega,mag_max,'r','LineWidth',1.5)
set(gca,'XScale','log')
set(gca,'YScale','log')

%%
close 4
magg = vpck(mag_max,omega);
figure(4)
wu = fitmag(magg);
[A,B,C,D] = unpck(wu);

Wu = zpk(tf(ss(A,B,C,D)))

%% Wu and Wt
close 5

figure(5)
hold on,grid on
bodemag(omega,Wu,'--b')
bodemag(omega,1/WT_inv,'b')


%% W1 and W2

W2 = (1+s/p3)^2/Tpo;
W2mod = tf(1/Tpo);

figure(2)
bodemag(omega,1/W2,'c')

W1 = 1/WS_inv1;
lambda = 1;
W1mod = zpk(minreal(W1*s^2/(s+lambda)^2,1e-3))

figure(1)
bodemag(omega,1/W1mod,'c')

%% h inf
[A,B,C,D] = linmod('generalized_plant')
M = ltisys(A,B,C,D)
M = sderiv(M,2,[1/(p3) 1])
M = sderiv(M,2,[1/(p3) 1])

[gopt,Gcmod1] = hinflmi(M,[1 1], 0, 0.01,[0 0 0])
[A,B,C,D] = ltiss(Gcmod1);
Gcmod = zpk(tf(ss(A,B,C,D)))

%% clearning
Gc = minreal(Gcmod*(s^2 + 1.904*s + 0.9727)*(1+s/5.945e06)/(s*(s+0.000308)),1e-3)
%% Ln
L = minreal(Gcmod*Ga*Gpn*Gs*Gf,1e-3)
Ln = minreal(Gcmod*Ga*Gpn*Gs*Gf,1e-3)

close 6
omega_L = logspace(-5,5,10000);

figure(6)
myngridst(Tpo,Spo)

hold on,grid on
nichols(L,omega_L,'b')

nichols(Ln,omega_L,'r')

Sn = minreal(1/(1+Ln),1e-3)
figure(1)
bodemag(omega,Sn,'r')

Tn = minreal(Ln/(1+Ln),1e-3)
figure(2)
bodemag(omega,Tn,'r')
%%  time simulation
close 7
t = linspace(0,0.2,2000);
yr = step(Tn*Kd,t);

figure(7)
plot(t,yr,'b','LineWidth',1.5)
yline(Kd*1.05)
yline(Kd*0.95)
xline(tr,'--')
xline(ts,'-')

yline(Kd)
%% da
close 8
da = Da0*ones(length(t),1);

yda = lsim(Gp*Sn,da,t);
figure(8)
hold on,grid on
plot(t,yda,'b','LineWidth',1.5)
yline(1.75e-2)
yline(-1.75e-2)

%% dp
t = linspace(0,1000,10000)
close 9
dp = Dp0*t;
ydp = lsim(Sn,dp,t);


figure(9)
hold on,grid on
plot(t,ydp,'b','LineWidth', 1.5)
yline(0)

%% u
close 10
t = linspace(0,0.5,5000);

u = step(Gc*Ga*Sn,t)

figure(10)
hold on,grid on
plot(t,u,'b','LineWidth',1.5)

%%
close 11
ds = as*sin(2000*t);

yds = lsim(Tn/Gs,ds,t);

figure(11)
hold on,grid on
plot(t,yds,'b','LineWidth',1.5)
yline(2e-4)
yline(-2e-4)

%% PR unstruct
close 13

figure(13)
bodemag(omega,Wu*Tn,'--b')
hold on,grid on
bodemag(omega,Sn/WS_inv,'--g')
bodemag(omega,Wu*Tn + Sn/WS_inv,'b')


%% Mu
close 14
WT_mu_RP4 = tf(1/M_HF);
[A,B,C,D] = linmod('N_scheme_RP4')
N = pck(A,B,C,D);

omega_RP4 = logspace(log10(2000),5,1000);

Nf = frsp(N,omega_RP4);

deltaset = [-1,1;-3,0;1,1];

mubounds = mu(Nf,deltaset)

figure(14)
vplot('liv,m',mubounds)

%% robust performance unstructured
close 15
WT_mu_RP4 = tf(1/M_HF);
[A,B,C,D] = linmod('N_scheme_RP4_U')
N = pck(A,B,C,D);

omega_RP4 = logspace(log10(2000),5,1000);

Nf = frsp(N,omega_RP4);

deltaset = [1,1;1,1];

mubounds = mu(Nf,deltaset)

figure(15)
vplot('liv,m',mubounds)

%%
K_min = 20;
K_max = 60;
K_n = (K_max + K_min)/2
WK = (K_max-K_min)/2

b_min = 48;
b_max = 120;
b_n = (b_max + b_min)/2
Wb = (b_max-b_min)/2
c = 0.707;
s = tf('s')
Gpn = K_n/(s*(1+2*c*s/b_n +(s/b_n)^2))
kp = dcgain(s*Gpn)


%%
close 16

[A,B,C,D] = linmod('N_scheme_RS')
N = pck(A,B,C,D);


Nf = frsp(N,omega);

deltaset = [-1,1;-3,0];

mubounds = mu(Nf,deltaset)

figure(16)
vplot('liv,m',mubounds)

%%
close 17

[A,B,C,D] = linmod('N_scheme_RS_U')
N = pck(A,B,C,D);

%%

Nf = frsp(N,omega);

deltaset = [1,1];

mubounds = mu(Nf,deltaset)

figure(17)
vplot('liv,m',mubounds)








