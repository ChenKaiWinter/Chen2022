function [d] = euclidean_dis(loc_1, loc_2)

d = 0;
for i=1:2
    d = d+(loc_1(i)-loc_2(i))^2;
end
d = sqrt(d);

end
