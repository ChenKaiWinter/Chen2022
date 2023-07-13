% Simulation 
%   "Beamforming Optimizaition for Full-Duplex NOMA Enabled Integrated Sensing 
%   and Communication" - Kai Chen et. al. 2022

% Kai Chen, 7/22/2020

%% Set up the system
P_t_dB = 20; % inital: 20, 25, 30, 35dBm transmit power at BS
P_u_dB = 10; % 10dBm transmit power at UUE

P_t = 10^(P_t_dB/10); % 100 mW
P_u = 10^(P_u_dB/10); % 10 mW

n_d_dB = -114; % -114dBm noise power sigma_d
n_u_dB = -114; % -114dBm noise power sigma_u
n_RSI_dB = -114; % inital: -114, -110, -106dBm noise power sigma_RSI

n_d = 10^(n_d_dB/10); 
n_u = 10^(n_u_dB/10);
n_RSI = 10^(n_RSI_dB/10);

n_0_dB = -104; % inital: -104dBm complex amplitude variance 
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

beta0_dB = -20; % inital: -20dBm
beta_0 = 10^(beta0_dB/10);

a_tx = zeros(N_t, length(angleSpace)); % transmit antennas steering vecter
for i = 1:N_t
    a_tx(i,:) = exp((1i * 2 * pi * TxAntLoc(i) / lambda) .* sin(angleSpace));
end
a_rx = zeros(N_r, length(angleSpace)); % receive antennas steering vecter
for i = 1:N_r
    a_rx(i,:) = exp((1i * 2 * pi * RxAntLoc(i) / lambda) .* sin(angleSpace));
end

%% location
tgt_num = 1; % number of Target
tgt_dis = [140];
tgt_agl = [0];
K_d = 2; % number of DUEs
DUE_dis = [10 10];
DUE_agl = [40 -40]; 
K_u = 4; % number of UUEs
UUE_dis = [60 80 100 120]; 
UUE_agl = [60 30 -30 -60];

a_T = zeros(N_t,tgt_num);
for k = 1:tgt_num
    a_T(:,k) = a_tx(:,90+tgt_agl(k));
end
a_R = zeros(N_r,tgt_num);
for k = 1:tgt_num
    a_R(:,k) = a_rx(:,90+tgt_agl(k));
end

% generate the user's coordinates
DUE_loc_x = zeros(1, K_d);
DUE_loc_y = zeros(1, K_d);
UUE_loc_x = zeros(1, K_u);
UUE_loc_y = zeros(1, K_u);
for i =1:K_d
    DUE_loc_x(1,i) = DUE_dis(i)*cos(DUE_agl(i));
end
for i =1:K_d
    DUE_loc_y(1,i) = DUE_dis(i)*sin(DUE_agl(i));
end
DUE_loc = [DUE_loc_x;DUE_loc_y];
for i =1:K_u
    UUE_loc_x(1,i) = UUE_dis(i)*cos(UUE_agl(i));
end
for i =1:K_u
    UUE_loc_y(1,i) = UUE_dis(i)*sin(UUE_agl(i));
end
UUE_loc = [UUE_loc_x;UUE_loc_y];

D_U_dis = zeros(K_u,K_d);
for i = 1:K_u
    for j = 1:K_d
        D_U_dis(i,j) = generate_dis(UUE_loc(:,i),DUE_loc(:,j));
    end
end

alpha = 3.8; % inital: 3.8
% beta_d = beta_0 * DUE_dis.^(-alpha);
% beta_u = beta_0 * UUE_dis.^(-alpha);
% beta_t = beta_0 * D_U_dis.^(-alpha);
beta_d = DUE_dis.^(-alpha);
beta_u = UUE_dis.^(-alpha);
beta_t = D_U_dis.^(-alpha);

%% begin
G_dB = 10; % inital: 10

[h_d, h_u, h_t] = generate_ch_Rci(G_dB, N_t, N_r, K_d, K_u, TxAntLoc, RxAntLoc, lambda, angleSpace, DUE_agl, UUE_agl, beta_d, beta_u, beta_t);
save('mat\ch.mat');

% load('mat\ch.mat');

load('mat\W_0.mat');

w = reshape(W_0, N_t.^2, 1); % inital transmit beamforming vecter full power

[~, A, B_i, R_ko_2] = generate_mt(a_R, a_T, n_RSI, n_u, n_0, N_t, N_r, w); % generate the matrix that we need

%% generate the R_{k_u,o} & MMSE receiver
R_ko_1 = zeros(N_r,N_r,K_u);
for i = 1:K_u
     if i == K_u
         R_ko_1(:,:,i) = zeros(N_r,N_r);break
     end
     for j = i+1:K_u
         R_ko_1(:,:,i) = R_ko_1(:,:,i) + P_u*h_u(:,j)*h_u(:,j)';
     end
end
R_ko = zeros(N_r,N_r,K_u);
for i = 1:K_u
    R_ko(:,:,i) = R_ko_1(:,:,i) + R_ko_2;
end
v_kc = zeros(N_r,1);
for i = 1:K_u
    v_kc(:,i) = MMSE(P_u, R_ko(:,:,i), h_u(:,i));
end

%% generate martix that I need
O = (n_0/(n_RSI+n_u))*(A*A');
Q_kd = zeros(N_t.^2,N_t.^2,K_d);
for i = 1:K_d
    Q_kd(:,:,i) = B_i(i)'*h_d(:,i)*h_d(:,i)'*B_i(i);
end
Q_i_1 = B_i(2)'*h_d(:,1)*h_d(:,1)'*B_i(2);
Q_i_2 = B_i(1)'*h_d(:,2)*h_d(:,2)'*B_i(1);
H_t = zeros(K_u,K_u,K_d);
for i = 1:K_d
    H_t(:,:,i) = P_u*h_t(:,i)*h_t(:,i)';
end

%% change this
gamma_d_dB = [14 14]; % 2 ~ 24dB inital: 14dB 
gamma_d = 10.^(gamma_d_dB/10);

%% CVX
cvx_begin SDP
variable R((N_t).^2, (N_t).^2) hermitian 
maximize (real(trace(O*R)))
subject to
real(trace(Q_kd(:,:,1)*R)) - gamma_d(1)*(real(trace(Q_i_1*R)) + real(trace(H_t(:,:,1))) + n_d) >= 0;
real(trace(Q_kd(:,:,2)*R)) - gamma_d(2)*(real(trace(Q_i_2*R)) + real(trace(H_t(:,:,2))) + n_d) >= 0;
trace(R) <= P_t;
R == semidefinite((N_t).^2);
cvx_end

if  rank(R)==1
    % eigenvalue decomposition: sqrt(dominant eigenvalue)*dominant eigenvector
    [V,D] = eigs(R);
    d = diag(D);
    [mv,mi]= max(d);
    w_dom = sqrt(d(mi))*V(:,mi);
    w_opt = w_dom;
else
%     % rank-one approximation for BQP
%     [V,D] = eigs(R);
%     d = abs(diag(D));
%     [mv,mi]= max(d);
%     w_r = sqrt(d(mi))*V(:,mi); % sqrt(dominant eigenvalue)*dominant eigenvector
%     w_opt = sign(w_r);  
    % Gaussian randomization
% Do randomization to get our feasible beamforming vectors w_cand
nRand = 100;
w_cand = zeros(N_t.^2, nRand);
W_cand = zeros(N_t.^2, N_t.^2, nRand);
clear uF;
% generate and scale to meet power constraints
for L = 1:nRand
    % generate nRand
    w_cand(:,L) = mvnrnd(zeros(N_t.^2,1),R) + 1i*mvnrnd(zeros(N_t.^2,1),R);
    % scale them so it adheres to norm constraint
    w_cand(:,L) = sqrt(trace(R))*w_cand(:,L)/sqrt(w_cand(:,L)'*w_cand(:,L));
    W_cand(:,:,L) = w_cand(:,L)*w_cand(:,L)';
end

index = 1;
for i = 1:nRand
if (real(trace(Q_kd(:,:,1)*W_cand(:,:,i)))-gamma_d(1)*(real(trace(Q_i_1*W_cand(:,:,i)))+real(trace(H_t(:,:,1)))+n_d)>=0)&...
        (real(trace(Q_kd(:,:,2)*W_cand(:,:,i)))-gamma_d(2)*(real(trace(Q_i_2*W_cand(:,:,i)))+real(trace(H_t(:,:,2)))+n_d)>=0)
    uF(:,1,index) = w_cand(:,i);
    index = index + 1;
end
end

numFeasible = size(uF,3);
W_opt = zeros(N_t.^2,N_t.^2,numFeasible); obj = zeros(1,numFeasible);
for i = 1:size(uF,3)
    W_opt(:,:,i) = (uF(:,1,i)*uF(:,1,i)');
    obj(i) = real(trace(O*W_opt(:,:,i)));
end
[maxVal,maxIndex] = max(obj);
w_opt = uF(:,1,maxIndex);
end

%% Before plot beampattern
% get the v_kc with SIC
[~, A, B_i, R_ko_2] = generate_mt(a_R, a_T, n_RSI, n_u, n_0, N_t, N_r, w_opt); % generate the matrix that we need
R_ko_1 = zeros(N_r,N_r,K_u);
for i = 1:K_u
     if i == K_u
         R_ko_1(:,:,i) = zeros(N_r,N_r);break
     end
     for j = i+1:K_u
         R_ko_1(:,:,i) = R_ko_1(:,:,i) + P_u*h_u(:,j)*h_u(:,j)';
     end
end
R_ko = zeros(N_r,N_r,K_u);
for i = 1:K_u
    R_ko(:,:,i) = R_ko_1(:,:,i) + R_ko_2;
end
% v_kc = zeros(N_r,1);
% for i = 1:K_u
%     v_kc(:,i) = MMSE(P_u, R_ko(:,:,i), h_u(:,i));
% end

W_opt = reshape(w_opt, N_t, N_t);

%% SINR_r
SINR_r = real(10*log10((n_0/n_RSI+n_u)*w_opt'*A'*A*w_opt));

%% SINR_u NOMA-c
% SINR_u with SIC
SINR_u = zeros(1,K_u);
for i = 1:K_u
SINR_u(:,i) = real(10*log10(P_u*h_u(:,i)'*inv(R_ko(:,:,i))*h_u(:,i)));
end
% SINR_u without SIC
SINR_u_c = zeros(1,K_u);
R_ko_c_sum = zeros(N_r,N_r);
R_ko_c = zeros(N_r,N_r,K_d);
for i = 1:K_u
    R_ko_c_sum = R_ko_c_sum + P_u*h_u(:,i)*h_u(:,i)';
end
for i = 1:K_u
    R_ko_c(:,:,i) = R_ko_c_sum - P_u*h_u(:,i)*h_u(:,i)' + R_ko_2;
end
for i = 1:K_u
SINR_u_c(:,i) = real(10*log10(P_u*h_u(:,i)'*inv(R_ko_c(:,:,i))*h_u(:,i)));
end

% SINR_d
SINR_d_1 = real(10*log10((abs(h_d(:,1)'*W_opt(:,1))).^2/((abs(h_d(:,1)'*W_opt(:,2))).^2 + P_u*h_t(:,1)'*h_t(:,1) + n_d)));
SINR_d_2 = real(10*log10((abs(h_d(:,2)'*W_opt(:,2))).^2/((abs(h_d(:,2)'*W_opt(:,1))).^2 + P_u*h_t(:,2)'*h_t(:,2) + n_d)));

%% Plot beampattern
TxBp = zeros(size(angleSpace));
for i = 1:length(angleSpace)
    TxBp(i) = abs(a_tx(:,i)' * (W_opt*W_opt') * a_tx(:,i));
end
figure; 
TxBp_l = plot(angleSpaceDeg, mag2db(TxBp), 'LineWidth', 1, 'Linestyle', '--'); 
tgt_l = line([tgt_agl, tgt_agl],[min(mag2db(TxBp)), max(mag2db(TxBp))],'Color', 'black', 'LineWidth', 0.5,'linestyle','--');
for i =1:K_d
DUE_l = line([DUE_agl(i),DUE_agl(i)],[min(mag2db(TxBp)), max(mag2db(TxBp))],'Color', 'magenta', 'LineWidth', 0.5,'linestyle','--');
end
hold on;
grid on;
xlabel('Angle Space [-90^\circ,90^\circ]'); ylabel('Magnitude (dB)')
title('ISAC Tranmit Beampattern'); 
legend([TxBp_l, tgt_l, DUE_l],'ISAC Tranmit Beampattern','Sensing Direction','Comm Direction'); 
axis([-90, 90, min(mag2db(TxBp)), max(mag2db(TxBp))]);