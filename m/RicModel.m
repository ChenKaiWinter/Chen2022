function H_ric = RicModel(G_dB, Nt_hor, Nt_ver, Nt, thetaCom, phiCom, lambda, x_t, y_t, z_t, K, beta)

%

G = 10^(G_dB/10);               % Rician factor


b_ric_vec = zeros(Nt,K);

for k = 1: K
    b_ric = zeros(Nt_hor , Nt_ver);
    for p = 1:Nt_hor  % x
        for q = 1:Nt_ver  % y
            position = [x_t(p) y_t z_t(q)]; %阵元的位置坐标
            direction_com = [sind(thetaCom)*cosd(phiCom(k)) 0 cosd(thetaCom)];
            b_ric(p,q) = exp(1i * 2 * pi * position * direction_com.'/ lambda);
        end
    end
    b_ric_vec(:,k) = reshape(b_ric,Nt,1);
end

H_ric = sqrt(beta) * (sqrt(G/(G+1)) * b_ric_vec + sqrt(1/(G+1)) * RayModel(Nt, K));
%H_ric = sqrt(G/(G+1)) * b_ric_vec + sqrt(1/(G+1)) * RayModel(Nt, K);
