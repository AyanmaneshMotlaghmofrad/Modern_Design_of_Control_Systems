clear all
close all
clc

%% Problem Definition
% Parameter uncertainties
K = [20, 60];
zeta_p = [0.65, 0.75];
w_np = [1.7, 2.5];
K_pn = mean(K);
zeta_pn = mean(zeta_p);
w_npn = mean(w_np);

% Plant 
s = tf('s');
G_pn = K_pn / (4.5 * (1 + 2*zeta_pn/w_npn*s + s^2/w_npn^2));
K_p = dcgain(G_pn);
G_s = 1;
G_a = -0.09;
D_a0 = 8.5e-3;
D_p0 = 3e-3;
a_p = 0; 
w_p = 0;
a_s = 1e-2;
w_s = 50;

%% Requirements translation
%S1
K_d = 1;
G_f = 1 / (K_d*G_s);

%S2
e_r = 3.5e-1;
R_0 = 1;
nu_1 = 1;
S_star_zero_1 = abs(e_r / (K_d*R_0)); % |S_star_zero| < 0.35

%S3
e_da = 1.75e-2;
nu_2 = 0;
S_star_zero_2 = abs(e_da / (K_p*D_a0)); % |S_star_zero| < 0.2316

%S4
e_dp = 1e-3;
nu_3 = 1;
S_star_zero_3 = abs(e_dp / D_p0); % |S_star_zero| < 0.3333

%S5
e_ds = 2e-4;
M_T_HF = e_ds * G_s / a_s;
M_T_HF_dB = 20*log10(M_T_HF); % M_T_HF_dB = -33.9794dB

%% Transient requirements
%S8
s_hat = 0.08;
zeta = abs(log(s_hat)) / (sqrt(pi^2 + log(s_hat)^2)); % zeta = 0.6266

%S6
t_r = 2.5;
w_n_tr = 1 / (t_r * sqrt(1 - zeta^2)) * (pi - acos(zeta)); % w_n > 1.1537rad/s 

%S7
t_s5 = 10;
w_n_ts5 = log(100/5) / (t_s5 * zeta); % w_n > 0.4781rad/s 

S_p = 2 * zeta * sqrt(2 + 4 * zeta^2 + 2 * sqrt(1 + 8 * zeta^2)) /...
    (sqrt(1 + 8 * zeta^2) + 4 * zeta^2 - 1);
T_p = 1 / (2 * zeta * sqrt(1 - zeta^2));

%% Choice of parameters
nu = max(nu_1, max(nu_2, nu_3)); % nu = 1
w_n_min = max(w_n_tr, w_n_ts5); % w_c_min = 1.1537rad/s

S = s * (s + 2*zeta*w_n_min) / (s^2 + 2*zeta*w_n_min*s + w_n_min^2);

%% Weighting funciton W_S_inv design
a = min(S_star_zero_1, S_star_zero_3);
S_star_zero = a*s; 

figure, hold on
bodemag(S, {1e-4, 1e+3}, '--g')
bodemag(S_star_zero, 'g') 
semilogx([1e-4, 1e+3], [20*log10(S_p), 20*log10(S_p)], 'r')

W_S_inv = a * s;
% bodemag(W_S_inv)

% Real zero
z1 = 0.01;
W_S_inv = W_S_inv * (1 + s / z1);
% bodemag(W_S_inv)

% Real pole
p1 = 0.033;
W_S_inv = W_S_inv * 1 / (1 + s/p1);
% bodemag(W_S_inv)

% Real zero and complex conjugate poles
z2 = 1;
p2 = sqrt(S_p*z1*z2/(a*p1));
zeta = 0.85;
W_S_inv = W_S_inv * (1 + s/z2) / (1 + 2*zeta*s/p2 + s^2/p2^2);
bodemag(W_S_inv, {1e-4, 1e+3})

grid on, hold off

%% Weighting funciton W_T_inv design
p1 = w_s / 10^((20*log10(T_p)-M_T_HF_dB)/40);
zeta = 0.8;
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
x_zeta_p = linspace(zeta_p(1), zeta_p(2), n);
x_w_np = linspace(w_np(1), w_np(2), n);
w = logspace(-3, 3, 1000);
mG_pnf = freqresp(G_pn, w);
count = 1;
mf = zeros(1000,1);

figure, 
for i = 1:n
    for j = 1:n
        for k = 1:n
            G_p = x_K(i) / (4.5 * (1 + 2*x_zeta_p(j)/x_w_np(k)*s + s^2/x_w_np(k)^2));
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
zeta = 0.6266;
w_c = w_n_min * sqrt(sqrt(1+4*zeta^4)-2*zeta^2);

% W_1
W_1 = 1 / W_S_inv;
figure, bodemag(W_1), grid on
W_1mod = zpk(minreal(W_1 * s / (s + 0.01*w_c)));
figure, hold on
bodemag(1 / W_1mod, 'b')
bodemag(S, {1e-4, 1e+3}, '--g')
bodemag(S_star_zero, 'g') 
semilogx([1e-4, 1e+3], [20*log10(S_p), 20*log10(S_p)], 'r')
grid on, hold off

% W_2
figure, bodemag(1 / W_T_inv), grid on, hold on
bodemag(W_u_m)
W_T = 1 / W_T_inv;
W_2 = W_T;
W_2mod = tf(1 / T_p);

mW_2mod = squeeze(freqresp(W_2mod, w));
mW_u = squeeze(freqresp(W_u_m, w));
mf = max(mW_2mod, mW_u);
magg = vpck(mf(:), w);

figure, 
W_2mod = fitmag(magg);
[A,B,C,D] = unpck(W_2mod);
[z,p,k] = ss2zp(A,B,C,D);
W_2mod = zpk(z,p,k);
figure, bodemag(1 / W_T_inv), grid on, hold on
bodemag(W_u_m), bodemag(W_2mod), grid on

%% 
[Am,Bm,Cm,Dm] = linmod('generalized_plant');
M = ltisys(Am,Bm,Cm,Dm);
p_WT = abs(pole(W_T_inv));
M = sderiv(M,2,[1/p_WT(1) 1]);
M = sderiv(M,2,[1/p_WT(2) 1]);
%M = sderiv(M,2,[1/p^2, 2*zeta/p, 1]);
[gopt,G_cmod] = hinflmi(M,[1,1],0,0.01,[0,0,0]);
[Ac,Bc,Cc,Dc] = ltiss(G_cmod);
G_cmod = ss(Ac,Bc,Cc,Dc);
G_cmod = minreal(zpk(G_cmod))

G_cmod = minreal(G_cmod * (s+0.008043) / s, 1e-3);
G_cmod = minreal(G_cmod * (1+2.703/9.315*s+s^2/9.315) / (1+2.404/9.122*s+s^2/9.122), 1e-3);
G_cmod = minreal(G_cmod * (1+s/1801), 1e-3)

%% Nichols Plot
L = G_cmod * G_a * G_pn * G_f * G_s;
S = 1 / (1+L);

figure,
myngridst(T_p, S_p), hold on
nichols(L)

%% Computation of errors for requirements
S_star = minreal(S * 1 / s);
T = L / (1 + L);
m_S = freqresp(S, w_p);
m_T = freqresp(T, w_s);

e_r = dcgain(S_star) * K_d
e_da = dcgain(S_star) * K_p * D_a0
e_dp = dcgain(S_star) * D_p0
e_ds = abs(m_T) * a_s / G_s

t_sim = 25;
out = sim("Simulator.slx");
figure,
plot(out.y.time, out.y.data)
yline(K_d+0.08*K_d, '--r', LineWidth=0.6)
yline(K_d+0.05*K_d, '--m', LineWidth=0.6)
yline(K_d-0.05*K_d, '--m', LineWidth=0.6)
xline(2.5, '--k', LineWidth=0.6)
xline(10, '--k', LineWidth=0.6)
grid on
