function chi = rt_kraus2chi(E)

s = size(E,1);
r = size(E,3);
e = zeros(s^2,r);

for k = 1:r
    e(:,k) = reshape(E(:,:,k), s^2, 1);
end

chi = e*e';

end

