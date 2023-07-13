function [h_d, h_u, h_t] = generate_ch_Ray( )
N_t = 8; % tranmit antennas number
N_r = 8; % receive antennas number
K_d = 2; % number of DUEs
DUE_dis = [50 50];
DUE_agl = [45 -45]; 
K_u = 4; % number of UUEs
UUE_dis = [20 40 60 80]; 
UUE_agl = [60 30 -30 -60];
%% position
for i=1:K_d
    x_DUE(i) = DUE_dis(i)*cos(90-DUE_agl(i));
end
for i=1:K_d
    y_DUE(i) = DUE_dis(i)*sin(90-DUE_agl(i));
end

for i=1:K_u
    x_UUE(i) = UUE_dis(i)*cos(90-UUE_agl(i));
end
for i=1:K_u
    y_UUE(i) = UUE_dis(i)*sin(90-UUE_agl(i));
end

DUE_loc = [x_DUE; y_DUE];
UUE_loc = [x_UUE; y_UUE];

%% Rayleigh channel
num = 100;

% downlink
h_d_all = zeros(N_t, K_d, num); % repeat 100 times
for k = 1:num
    for i = 1:K_d
        for j = 1:N_t
            h_d_all =  ( randn(j, i, k) + 1i*randn(j, i, k) )/sqrt(2);
        end
    end
end
h_d_s = zeros(N_t, K_d);
for i = 1:num
    h_d_s = h_d_s + h_d_all(:,:,i);
end
h_d_s = h_d_s/num;
for i = 1:K_d
    h_d(:,i) = sqrt(DUE_dis(i).^(-3.8)).*h_d_s(:,i);
end

% uplink
h_u_all = zeros(N_t, K_u, num);
for k = 1:num
    for i = 1:K_u
        for j = 1:N_t
            h_u_all =  ( randn(j, i, k) + 1i*randn(j, i, k) )/sqrt(2);
        end
    end
end
h_u_s = zeros(N_r,K_u);
for i = 1:num
    h_u_s = h_u_s + h_u_all(:,:,i);
end
h_u_s = h_u_s/num;
for i = 1:K_u
    h_u(:,i) = sqrt(UUE_dis(i).^(-3.8)).*h_u_s(:,i);
end

% DUE-UUE
for i = 1:K_d
    for j = 1:K_u
        d_t(j,i) = euclidean_dis(DUE_loc(:,i), UUE_loc(:,j));
    end
end

save('mat\d','d_t')

h_t_all = zeros(K_u, K_d, num); % repeat 100 times
for k = 1:num
    for i = 1:K_d
        for j = 1:K_u
            h_t_all = ( randn(j, i, k) + 1i*randn(j, i, k) )/sqrt(2);
        end
    end
end
h_t_s = zeros(K_u, K_d);
for i = 1:num
    h_t_s = h_t_s + h_t_all(:,:,i);
end
h_t_s = h_t_s/num;
for i = 1:K_d
    for j = 1:K_u
        h_t(j,i) = sqrt(d_t(j,i).^(-3.8)).*h_t_s(j,i);
    end
end

save('mat\ch.mat', 'h_d', 'h_u', 'h_t')

end