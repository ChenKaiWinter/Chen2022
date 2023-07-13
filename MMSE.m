function [vk_c] = MMSE(P_u, R_ko, h_u)

vk_c = sqrt(P_u)*inv(R_ko)*h_u;  

end