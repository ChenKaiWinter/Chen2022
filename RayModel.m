function H_ray = RayModel(Nt,K)

for i = 1:K
H_ray(:,i) = (randn(Nt,1)+1i*randn(Nt,1))/sqrt(2);
end