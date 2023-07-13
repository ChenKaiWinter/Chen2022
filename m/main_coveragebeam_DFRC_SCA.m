% ---------------------------------------------
% Main function of coverage-aware beamforming optimization
% The antenna arrays are located at X-Z plate
% 
% ---------------------------------------------

clc
clear all
%% URA parameter setup

Nt_hor = 8;
Nt_ver = 8;

Nr_hor = 4;                    
Nr_ver = 4;

Nt = Nt_hor * Nt_ver;            % number of Tx
Nr = Nr_hor * Nr_ver;            % number of Rx

fc = 5e9;           % center freq 5GHz
c = 3e8;            % speed of light
lambda = c/fc;      % wavelength   
spacing = lambda/2; % antenna spacing
B = 1e6;            % 1MHz
beta0_dB = -20;
T = 10^-3;

beta_s_dB = -40;
beta_s = 10^(beta_s_dB/10);
beta_s2 = beta_s^2;


% resolution 0.5
thetaSpace = linspace(0,180,361);  % range of elevation angle [0, pi]
phiSpace = linspace(0,180,361);  % range of azimuth angle     [0, pi]

x_t = spacing*[0:Nt_hor-1]; % location of Tx
y_t = 0;
z_t = spacing*[0:Nt_ver-1];

x_r = spacing*[0:Nr_hor-1];  % location of Rx
y_r = 0; 
z_r = spacing*[0:Nr_ver-1];

K = 2;          % number of UEs
thetaCom = 90;  % location of communication user
phiCom = [20,160];
beta_0 = 10^(beta0_dB/10); % pathloss
d = 30; % distance between BS-1 and UEs
alpha = 3.2;

beta = beta_0 * d^(-alpha);


thetaTar = 90;  % location of interested region
phiTar = [45,135];  % The minimum and maximum of azimuth angle of interested region
phiTarUnif = linspace(45,135,91);   % Uniformly saperate the angle


Pt_dBm = 20;
Pt = 10^(Pt_dBm/10);
%Pt = 1e-3 * Pt_mW;   %Pt mWatt to Watt

G_dB = 10;                  % Rician factor dB

h = RicModel(G_dB, Nt_hor, Nt_ver, Nt,thetaCom, phiCom, lambda, x_t, y_t, z_t, K, beta); % Communication channel Rician channel

gamma_th_dB = 10; % thershold of SINR 10dB
gamma_th = 10^(gamma_th_dB/10);


sigma2_dBm = -110;              % AWGN
sigma2 = 10^(sigma2_dBm/10); % unit Watt
sigma = sqrt(sigma2);

%% Sensing Region Settings %%

% z = 0 ---------------------------------
o = [0, 0];              %coordinate of BS-1
o_prime = [0, 200];      %coordinate of BS-2

% we set the endpoints of the rectangle (interested region) are located at
% (60*sqrt(3), 60, 0), (30*sqrt(3), 60, 0), (30*sqrt(3), 90, 0),
% (60*sqrt(3), 90, 0)----------------------------------

N = 50;  % divide the region into 20*20 parts 

x_q = linspace(-25, 25, N);                      % width 50 * length 60  (-25, 25 * sqrt(3)), (25, 25 * sqrt(3)), (-25, 60 + 25 * sqrt(3))
y_q = linspace(25, 75, N);

r1 = zeros(length(x_q),length(y_q));
r2 = zeros(length(x_q),length(y_q));
phiTarDeg = zeros(length(x_q),length(y_q));

for i = 1:length(x_q)
    for j = 1:length(y_q)
        r1(i,j) = norm([x_q(i),y_q(j)]).^2;
        r2(i,j) = norm([x_q(i),y_q(j)]-o_prime).^2;
        if x_q(i) >= 0 
            phiTarDeg(i,j) = acotd(x_q(i)/y_q(j)); 
        else
            phiTarDeg(i,j) = acotd(x_q(i)/y_q(j)) + 180;    % angle of each coordinate in the region
        end
    end
end
r = r1 .* r2; % product of distance d_1^2 * d_2^2
eta = r ./ Nr;

phiTarDeg = reshape(phiTarDeg, N*N, 1);
eta = reshape(eta, N*N, 1);

paraReg_temp1 = cell(N*N,2); % cell for saving coordinates of 

paraReg_temp1(:,1) = mat2cell(phiTarDeg, ones(1,N*N), 1);
paraReg_temp1(:,2) = mat2cell(eta, ones(1,N*N), 1);



paraReg_temp2 = cell2mat(paraReg_temp1(:,1));
[~,ind] = sort(paraReg_temp2);
paraReg = paraReg_temp1(ind,:);

% steering vector for each direction b(q_l)

%phiTarReg = reshape (phiTarDeg_temp, N * N, 1);
phiTarDeg = cell2mat(paraReg(:,1));
b_tar = zeros(Nt, length(phiTarDeg));     % steering vector of targeted

for i = 1 :length(phiTarDeg)
    b_tar_temp = zeros(Nt_hor , Nt_ver);
    for p=1:Nt_hor  % x
        for q=1:Nt_ver  % y
            position=[x_t(p) y_t z_t(q)]; %阵元的位置坐标
            direction_tar=[sind(thetaTar)*cosd(phiTarDeg(i)) sind(thetaTar)*sind(phiTarDeg(i)) cosd(thetaTar)];
            b_tar_temp(p,q) = exp(1i * 2 * pi * position * direction_tar.'/ lambda);
        end
    end
    b_tar(:,i) = reshape(b_tar_temp, Nt, 1);
end
% steering vector of targeted region (angle) ------------
b_tarAngle = zeros(Nt, length(phiTarUnif));     

for i = 1 :length(phiTarUnif)
    b_tar_temp = zeros(Nt_hor , Nt_ver);
    for p=1:Nt_hor  % x
        for q=1:Nt_ver  % y
            position=[x_t(p) y_t z_t(q)]; %阵元的位置坐标
            direction_tar=[sind(thetaTar)*cosd(phiTarUnif(i)) sind(thetaTar)*sind(phiTarUnif(i)) cosd(thetaTar)];
            b_tar_temp(p,q) = exp(1i * 2 * pi * position * direction_tar.'/ lambda);
        end
    end
    b_tarAngle(:,i) = reshape(b_tar_temp, Nt, 1);
end

% steering vector of all directions -------------
b_All = zeros(Nt, length(phiSpace));     

for i = 1 :length(phiSpace)
    b_tar_temp = zeros(Nt_hor , Nt_ver);
    for p=1:Nt_hor  % x
        for q=1:Nt_ver  % y
            position=[x_t(p) y_t z_t(q)]; %阵元的位置坐标
            direction_tar=[sind(thetaTar)*cosd(phiSpace(i)) sind(thetaTar)*sind(phiSpace(i)) cosd(thetaTar)];
            b_tar_temp(p,q) = exp(1i * 2 * pi * position * direction_tar.'/ lambda);
        end
    end
    b_All(:,i) = reshape(b_tar_temp, Nt, 1);
end


%% Beamforming Optimization %%
[W_opt1,zeta_opt] = SCA_sensing(Nt, h, Pt, b_tar, paraReg, K, sigma, gamma_th); % point discretization
W_opt2 = SCA_sensingAngle(Nt, h, Pt, b_tarAngle, K, sigma, gamma_th, eta); % point discretization

%% Plot beampattern

% beampattern -------------------------
beamP1 = zeros(length(phiSpace),1);
beamP2 = zeros(length(phiSpace),1);


for i = 1:length(phiSpace)
    beamP1(i) = abs(ctranspose(b_All(:,i)) * W_opt1 * ctranspose(W_opt1) * b_All(:,i))/real(trace(W_opt1 * ctranspose(W_opt1)));
    beamP2(i) = abs(ctranspose(b_All(:,i)) * W_opt2 * ctranspose(W_opt2) * b_All(:,i))/real(trace(W_opt2 * ctranspose(W_opt2)));
end

% SNR_zeta --------------------------------

SNR = zeros(length(zeta_opt),1);
SNR_dB = zeros(length(zeta_opt),1);

for i = 1:length(zeta_opt)
    SNR(i) = zeta_opt(i) * T * B * beta_s2/sigma2;
    SNR_dB(i) = 10 * log10(SNR(i));
end


%% plot figure

figure(1)
plot(phiSpace, mag2db(beamP1),'LineWidth', 2);
hold on
plot(phiSpace, mag2db(beamP2),'LineWidth', 2);
hold on

figure(2)
plot(SNR_dB,'LineWidth', 2);


%p_c = line([phiSpace(61) phiSpace(61)], [min(mag2db(beamP1)) max(mag2db(beamP1))], 'Color', 'magenta', 'LineStyle', '--');
p_c = line([phiSpace(91) phiSpace(91)], [min(mag2db(beamP1)) max(mag2db(beamP1))], 'Color', 'magenta', 'LineStyle', '--');






