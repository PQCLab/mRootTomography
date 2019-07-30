function d = rt_dm_theory(rho, M, n)
%rt_DM_THEORY Summary of this function goes here
%   Detailed explanation goes here

% Chi purification
r = rank(rho);
[U,D] = svd(rho);
U = U(:,1:r);
D = D(1:r,1:r);
c = U*sqrt(D);

% Information matrix
m = length(n);
s = size(rho,1);
H = 0;
for j = 1:m
    Hj = 0;
    for i = 1:s
        p = real(trace(M{j}(:,:,i)*rho))+1e-8;
        a = M{j}(:,:,i)*c;
        a = a(:);
        a = [real(a); imag(a)]; 
        Hj = Hj + (a*a')/p;
    end
    H = H + n(j)*Hj;
end
H = 2*H;

% Finding d
h = svd(H);
h(h < 1) = [];
d = 1./(2*h);

end

