function [R_p, A, B_i, R_ko_2] = generate_mt(a_R, a_T, n_RSI, n_u, n_0, N_t, N_r, w)

% R_p
R_p = n_RSI*eye(N_r*N_t)+n_u*eye(N_r*N_t);

% A(\theta_0)
A = kron(eye(N_t), a_R*a_T.');

% B_i
B_i = @(i) [zeros(N_t, N_t*(i-1)), eye(N_t), zeros(N_t, N_t*(N_t-i))];

% B
B = zeros(N_t, (N_t).^2);
for i = 1:N_t
    B =  B + B_i(i); % B
end
% C(\theta_0)
C = a_R*a_T.'*B;
% second parts of the Rk_o
R_ko_2 = n_0*C*w*w'*C' + n_RSI*eye(N_r) + n_u*eye(N_r);

end