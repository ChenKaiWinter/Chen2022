%% generate inital W_0
P_t_dB = 20; % 20dBm transmit power at BS
P_u_dB = 17; % 17dBm transmit power at UUE

P_t = 10^(P_t_dB/10); % 100 mW
P_u = 10^(P_u_dB/10); % 50 mW

n_d_dB = -114; % -114dBm noise power sigma_d
n_u_dB = -114; % -114dBm noise power sigma_u
n_RSI_dB = -114; % -114dBm noise power sigma_RSI

n_d = 10^(n_d_dB/10); 
n_u = 10^(n_u_dB/10);
n_RSI = 10^(n_RSI_dB/10);

n_0_dB = -94; % -84dBm complex amplitude variance 
n_0 = 10^(n_0_dB/10);

N_t = 8; % tranmit antennas number
N_r = 8; % receive antennas number

f_c = 5e9;
c = 3e8;
lambda = c/f_c;
spacing = lambda/2; % half wavelength
TxAntLoc = spacing*[0:N_t-1]; % transmit antennas
RxAntLoc = spacing*[0:N_r-1]; % receive antennas

M = 181; % number of points on the angular grid
angleSpace = linspace(-pi/2, pi/2, M);
angleSpaceDeg = linspace(-90, 90, M);

beta0_dB = -20;
beta_0 = 10^(beta0_dB/10);

a_tx = zeros(N_t, length(angleSpace)); % transmit antennas steering vecter
for i = 1:N_t
    a_tx(i,:) = exp((1i * 2 * pi * TxAntLoc(i) / lambda) .* sin(angleSpace));
end
a_rx = zeros(N_r, length(angleSpace)); % receive antennas steering vecter
for i = 1:N_r
    a_rx(i,:) = exp((1i * 2 * pi * RxAntLoc(i) / lambda) .* sin(angleSpace));
end

tgt_num = 1; % number of Target
tgt_dis = [50];
tgt_agl = [0];
K_d = 2; % number of DUEs
DUE_dis = [10 10];
DUE_agl = [45 -45]; 
K_u = 4; % number of UUEs
UUE_dis = [100 120 140 160]; 
UUE_agl = [60 30 -30 -60];

a_T = zeros(N_t,tgt_num);
for k = 1:tgt_num
    a_T(:,k) = a_tx(:,90+tgt_agl(k));
end
a_R = zeros(N_r,tgt_num);
for k = 1:tgt_num
    a_R(:,k) = a_rx(:,90+tgt_agl(k));
end

alpha = 3.8;
beta_d = beta_0 * DUE_dis.^(-alpha);
beta_u = beta_0 * UUE_dis.^(-alpha);

gamma_d_dB = [2 2];
gamma_d = 10.^(gamma_d_dB/10);

G_dB = 20;

[h_d, h_u] = generate_ch_Rci(G_dB, N_t, N_r, K_d, K_u, TxAntLoc, RxAntLoc, lambda, angleSpace, DUE_agl, UUE_agl, beta_d, beta_u);

while true
w = sqrt(P_t/(2*N_t.^2))*randn((N_t).^2, 1)+sqrt(P_t/(2*N_t.^2))*1i*randn((N_t).^2, 1);
W = reshape(w,N_t,N_t);
if (trace(W*W')<=100) && (10*log10((abs(h_d(:,1)'*W(:,1))).^2/((abs(h_d(:,1)'*W(:,2))).^2  + n_d))>=gamma_d(1)) && (10*log10((abs(h_d(:,2)'*W(:,2))).^2/((abs(h_d(:,2)'*W(:,1))).^2  + n_d))>=gamma_d(2))
    W_0 = W;break
end
end
save('mat\W_0.mat','W_0')
% I know SINR_d is incorrect, but I just don't want to do this anymore. 