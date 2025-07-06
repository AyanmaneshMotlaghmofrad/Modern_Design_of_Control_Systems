clc
clear variables
close all
format compact

%% defining the system
s = tf('s');
Gp = 100/(s^3+5.5*s^2+4.5*s);
p = 1;
nu = 1;
s_star_0 = 64.2;
wp = 0.03;
ws = 60;
wn = 1.42;
Tp0 = 1.095;
Sp0 = 1.41;
zeta = 0.54;

%% Ws_meno_1
close all
S2 = s*(s+2*zeta*wn)/(s^2+2*zeta*wn*s+wn^2);
S_0 = s^(p+nu)*s_star_0;
omega = logspace(-4,2,1000);

% plots
figure(1)
hold on,
bodemag(S2,omega,'--b')
bodemag(S_0,{0.0001,1},'--r')
bodemag(tf(Sp0),omega,'--g')
bodemag(tf(10^(-30/20)),{0.0001,wp},'--m')
grid on

p1 = 0.011;
z1 = 0.65;
p2 = sqrt((Sp0*z1)/(p1*s_star_0));
Ws_meno_1 = s^(p+nu)*s_star_0*((1+s/z1))/((1+s/p1)*(1+(2*0.75*s/p2)+(s/p2)^2));

bodemag(Ws_meno_1,omega,'-k')

%% close all
figure(2)

T2 = 1/(1+(2*zeta*s/wn)+(s/wn)^2);
hold on,
bodemag(T2,omega,'--b')
bodemag(tf(Tp0),omega,'--g')
bodemag(tf(10^(-48/20)),{ws,100},'--m')
grid on
p3 = ws*10^((-48-20*log10(Tp0))/40)
T = Tp0/(1+(1.414*s/p3)+(s/p3)^2);
bodemag(T,omega,'g')
