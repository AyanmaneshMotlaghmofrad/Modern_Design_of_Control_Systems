clear all
close all
clc

%% Problem Definition
% Parameter uncertainties
K = [9, 16];
pole1 = [0.55, 1.05];
pole2 = [1.9, 3.1];
K_pn = mean(K);
p1_n = mean(pole1);
p2_n = mean(pole2);

% Plant 
s = tf('s');
G_pn = K_pn / (s * (1 + s/p1_n) * (1 + s/p2_n));
K_p = dcgain(s*G_pn);
G_s = 2;
G_a = 0.38;
D_a0 = 5.5e-3;
D_p0 = 0;
a_p = 2e-2;
w_p = 0.02;
a_s = 1e-1;
w_s = 40;

%% Requirements translation
%S1
K_d = 4;
G_f = 1 / (K_d*G_s);

%S2
e_r = 1.5e-1;
R_0 = 1;
nu_1 = 0;
S_star_zero_1 = e_r / (K_d*R_0); % |S_star_zero| < 0.0375

%S3
e_da = 5.8;
nu_2 = 1;
S_star_zero_2 = e_da / (K_pn*D_a0); % |S_star_zero| < 84.3636

%S4
e_dp = 3.6e-4;
M_S_LF = e_dp / a_p;
M_S_LF_dB = 20*log10(M_S_LF); % M_S_LF_dB = -34.8945
w_l = w_p * 10^(-M_S_LF_dB/40);
w_c_1 = 2*w_l; % w_c > 0.2981

%S5
e_ds = 1.25e-4;
M_T_HF = e_ds * G_s / a_s;
M_T_HF_dB = 20*log10(M_T_HF); % M_T_HF_dB = -52.0412

%% Transient requirements
%S8
s_hat = 0.12;
zeta = abs(log(s_hat)) / (sqrt(pi^2 + log(s_hat)^2)); % zeta = 0.5594

%S6
t_r = 2.5;
w_n_tr = 1 / (t_r * sqrt(1 - zeta^2)) * (pi - acos(zeta)); % w_n > 1.0445rad/s

%S7
t_s5 = 5;
w_n_ts5 = log(100/5) / (t_s5 * zeta); % w_n > 1.0710rad/s

S_p = 2 * zeta * sqrt(2 + 4 * zeta^2 + 2 * sqrt(1 + 8 * zeta^2)) /...
    (sqrt(1 + 8 * zeta^2) + 4 * zeta^2 - 1);
T_p = 1 / (2 * zeta * sqrt(1 - zeta^2));

%% Choice of parameters
nu = max(nu_1, nu_2); % nu = 1
w_n_min = max(w_n_tr, w_n_ts5); % w_c_min = 1.0710rad/s

S = s * (s + 2*zeta*w_n_min) / (s^2 + 2*zeta*w_n_min*s + w_n_min^2);

%% Weighting funciton W_S_inv design
a = S_star_zero_2;
S_star_zero = a*s^2; 

figure, hold on
bodemag(S, {1e-4, 1e+3}, '--g')
bodemag(S_star_zero, 'g') 
semilogx([1e-4, w_p], [M_S_LF_dB, M_S_LF_dB], 'r', LineWidth=1.5)
semilogx([w_p, w_p], [M_S_LF_dB, 0], 'r', LineWidth=1.5)
semilogx([1e-4, 1e+3], [20*log10(S_p), 20*log10(S_p)], 'r')

W_S_inv = a * s^2;
% bodemag(W_S_inv)

% Real pole
p1 = 0.0125;
W_S_inv = W_S_inv * 1 / (1 + s/p1);
% bodemag(W_S_inv)

% Real zero and complex conjugate poles
z2 = 0.35;
p2 = sqrt(S_p*z2/(a*p1));
zeta = 0.8;
W_S_inv = W_S_inv * (1 + s/z2) / (1 + 2*zeta*s/p2 + s^2/p2^2);
bodemag(W_S_inv, {1e-4, 1e+3})

grid on, hold off

%% Weighting funciton W_T_inv design
p1 = w_s / 10^((20*log10(T_p)-M_T_HF_dB)/40);
zeta = 0.72;
W_T_inv = T_p / (1 + 2*zeta*s/p1 + s^2/p1^2);

figure, hold on
bodemag(W_T_inv, {1e-3, 1e+3}, 'g')
semilogx([w_s, 1e+3], [M_T_HF_dB, M_T_HF_dB], 'r', LineWidth=1.5)
semilogx([w_s, w_s], [M_T_HF_dB, 0], 'r', LineWidth=1.5)
semilogx([1e-3, 1e+3], [20*log10(T_p), 20*log10(T_p)], 'r')

grid on, hold off

%% Weighting function W_u
n = 10;
x_K = linspace(K(1), K(2), n);
x_p1 = linspace(pole1(1), pole1(2), n);
x_p2 = linspace(pole2(1), pole2(2), n);
w = logspace(-3, 3, 1000);
mG_pnf = freqresp(G_pn, w);
count = 1;
mf = zeros(1000,1);

figure, 
for i = 1:n
    for j = 1:n
        for k = 1:n
            G_p = x_K(i) / (s * (1 + s/x_p1(j)) * (1 + s/x_p2(k)));
            mG_pf = freqresp(G_p, w);
            Delta_mf = abs((mG_pf ./ mG_pnf) - 1);
            mag = squeeze(Delta_mf);
            loglog(w, mag), hold on
            mf = max(mf, mag);
        end
    end
end

loglog(w, mf, '--c', 'LineWidth', 1.5)
magg = vpck(mf(:), w);

figure,
W_u_m = fitmag(magg);
[A,B,C,D] = unpck(W_u_m);
[z,p,k] = ss2zp(A,B,C,D);
W_u_m = zpk(z,p,k);

%% Robust Control
zeta = 0.5594;
w_c = w_n_min * sqrt(sqrt(1+4*zeta^4)-2*zeta^2);

% W1
W_1 = zpk(1 / W_S_inv); 
W_1mod = zpk(minreal(W_1 * s^2 / (s + 0.001*w_c)^2));
figure, hold on
bodemag(1 / W_1mod, 'b')
bodemag(S, {1e-4, 1e+3}, '--g')
bodemag(S_star_zero, 'g') 
semilogx([1e-4, w_p], [M_S_LF_dB, M_S_LF_dB], 'r', LineWidth=1.5)
semilogx([w_p, w_p], [M_S_LF_dB, 0], 'r', LineWidth=1.5)
semilogx([1e-4, 1e+3], [20*log10(S_p), 20*log10(S_p)], 'r')
grid on, hold off

% W2
W_T = zpk(1 / W_T_inv);
figure, hold on
bodemag(W_T), bodemag(W_u_m), grid on % No overlap 
W_2 = W_T;
W_2mod = tf(1 / T_p);

%%
[Am,Bm,Cm,Dm] = linmod('generalized_plant');
M = ltisys(Am,Bm,Cm,Dm);
p_WT = abs(pole(W_T_inv));
M = sderiv(M,2,[1/p_WT(1) 1]);
M = sderiv(M,2,[1/p_WT(2) 1]);
[gopt,G_cmod] = hinflmi(M,[1,1],0,0.01);
[Ac,Bc,Cc,Dc] = ltiss(G_cmod);
G_cmod = ss(Ac,Bc,Cc,Dc);
G_cmod = minreal(zpk(G_cmod))

% Controller Cancellations
G_cmod = minreal(G_cmod * (s+0.000853)/s, 1e-3);
G_cmod = minreal(G_cmod * (s+0.0007006) / (s+1.082e-05), 1e-3);
G_cmod = minreal(G_cmod * (1 + s/2.887e04), 1e-3)

%% Nichols Plot
L = G_cmod * G_a * G_pn * G_f * G_s;
S = 1 / (1+L);

figure,
myngridst(T_p, S_p), hold on
nichols(L)

%% Computation of errors for requirements
S_star = minreal(S * 1 / s^2);
T = L / (1 + L);
m_S = freqresp(S, w_p);
m_T = freqresp(T, w_s);

e_r = 0 % because nu=1
e_da = dcgain(S_star) * K_p * D_a0
e_dp = abs(m_S) * a_p
e_ds = abs(m_T) * a_s / G_s

t_sim = 25;
out = sim("Simulator.slx");
figure,
plot(out.y.time, out.y.data)
yline(K_d+0.12*K_d, '--r', LineWidth=0.6)
yline(K_d+0.05*K_d, '--m', LineWidth=0.6)
yline(K_d-0.05*K_d, '--m', LineWidth=0.6)
xline(2.5, '--k', LineWidth=0.6)
xline(5, '--k', LineWidth=0.6)
grid on

