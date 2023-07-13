function [H_ric_d, H_ric_u, H_t] = generate_ch_Rci(G_dB, N_t, N_r, K_d, K_u, TxAntLoc, RxAntLoc, lambda, angleSpace, DUE_agl, UUE_agl, beta_d, beta_u, beta_t)

G = 10^(G_dB/10); % Rician factor

tx_ric_vec = zeros(N_t,K_d);
rx_ric_vec = zeros(N_t,K_u);

for k = 1: K_d
    tx_ric = zeros(N_t,length(angleSpace));
    for i = 1: N_t
        tx_ric(i,:) = exp((1i * 2 * pi * TxAntLoc(i) / lambda) .* sin(angleSpace));
    end
end
for k = 1: K_u
    rx_ric = zeros(N_t,length(angleSpace));
    for i = 1: N_r
        rx_ric(i,:) = exp((1i * 2 * pi * RxAntLoc(i) / lambda) .* sin(angleSpace));
    end
end

for k = 1:K_d
    tx_ric_vec(:,k) = tx_ric(:,90+DUE_agl(k));
end

for k = 1:K_u
    rx_ric_vec(:,k) = rx_ric(:,90+UUE_agl(k));
end

H_ray_d = zeros(N_t,K_d);
for i = 1:K_d
H_ray_d(:,i) = (randn(N_t,1)+1i*randn(N_t,1))/sqrt(2);
end
H_ray_u = zeros(N_r,K_u);
for i = 1:K_u
H_ray_u(:,i) = (randn(N_r,1)+1i*randn(N_r,1))/sqrt(2);
end
H_ray_t = zeros(K_u,K_d);
for i = 1:K_d
H_ray_t(:,i) = (randn(K_u,1)+1i*randn(K_u,1))/sqrt(2);
end

H_ric_d = zeros(N_t,K_d);
for i = 1 : K_d
% H_ric_d(:,i) = sqrt(beta_d(i)) * (sqrt(G/(G+1)) * tx_ric_vec(:,i) + sqrt(1/(G+1)) * H_ray_d(:,i));
H_ric_d(:,i) = tx_ric_vec(:,i);
end
H_ric_u = zeros(N_t,K_u);
for i = 1 : K_u
% H_ric_u(:,i) = sqrt(beta_u(i)) * (sqrt(G/(G+1)) * rx_ric_vec(:,i) + sqrt(1/(G+1)) * H_ray_u(:,i));
H_ric_u(:,i) = rx_ric_vec(:,i);
end
% For downlink users and uplink users, 
% the angle information is no longer important. 
% Here, only the Rayleigh channel with path loss is considered.
% So what?
H_t = zeros(K_u,K_d);
for i = 1 : K_u
    for j = 1 : K_d
        H_t(i,j) = sqrt(beta_t(i,j)) * H_ray_t(i,j);
    end
end
