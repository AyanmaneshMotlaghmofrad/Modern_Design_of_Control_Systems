clc
clear variables
close all
format compact

%%
s = tf('s');
kmin = 9;
kmax = 16;
kn = (kmin+kmax)/2;

p1min = 0.55;
p1max = 1.05;
p1n = (p1min + p1max)/2;

p2min = 1.9;
p2max = 3.1;
p2n = (p2min + p2max)/2;

Gpn = kn/(s*(1 + s/p1n) * (1 + s/p2n));

Tp_max = 1.05;
Sp_max = 1.36;
wc_min = 0.66;%non serve
wc_max = 1.4;%non serve
S_star0 = 0.15;
nu = 0;
p = 1;
M_LF = -32;
 
%%
close 1
wn2 = 0.91;
zeta2 = 0.59;

S2 = s*(s + 2*zeta2*wn2)/(s^2 + 2*zeta2*wn2*s + wn2^2);
S_0 = s*S_star0;

figure(1)
grid on
bodemag(S2,{0.0001,1000},'--b')
hold on
bodemag(S_0,{0.0001,10},'--g')
bodemag(tf(Sp_max),{0.0001,1000},'--r')
bodemag(tf(10^(M_LF/20)),{0.0001,0.02})

z1 = 0.003;
p1 = 0.03;
z2 = 0.8;
p2 = sqrt(z1*z2*Sp_max/(S_star0*p1));
zeta =0.85;
W_S_inv = S_0 * (1 + s/z1)*(1 + s/z2)/((1+(2*zeta*s/p2)+(s/p2)^2)*(1+s/p1));
bodemag(W_S_inv,{0.0001,1000},'-k')

%% Weighting function on T: W_T
close 2
T2 = wn2^2/(s^2 + 2*zeta2*wn2*s + wn2^2);
p3 = 6.18;

figure(2)
bodemag(T2,{0.001,1000},'--b')
hold on, grid on
bodemag(tf(Tp_max),{0.001,1000},'--r')
bodemag(tf(10^(-32/20)),{40,1000},'--k')

W_T_inv = Tp_max/(1+(2*0.707*s/p3)+(s/p3)^2);
bodemag(W_T_inv,{0.001,1000},'k')

%% multiplicative uncertainty

close all

omega = logspace(-4,+4,1000);
Gpn_mag = squeeze(freqresp(Gpn,omega));

mag_max = -inf * ones(length(omega),1);
figure(3)
hold on
for k = linspace(9,16,10)
    for p1 = linspace(0.55,1.05,10)
        for p2 = linspace(1.9,3.1,10)
            Gp = k/(s*(1 + s/p1) * (1 + s/p2));
            Gp_mag = squeeze(freqresp(Gp,omega));
            delta = abs(Gp_mag./Gpn_mag -1);
            
            plot(omega,delta)
            mag_max = max(mag_max,delta);
        end
    end
end

plot(omega, mag_max, 'LineWidth', 2, 'Color', 'r');
set(gca,'XScale','log')
set(gca,'YScale','log')

%% Wu
figure(4)
magg = vpck(mag_max,omega);
Wu = fitmag(magg);
[A,B,C,D] = unpck(Wu);

%%
sys = ss(A,B,C,D);
Wu = tf(sys);

Wu_mag = abs(squeeze(freqresp(Wu,omega)));

figure(3)
plot(omega,Wu_mag,'g','LineWidth',1.4)

%%