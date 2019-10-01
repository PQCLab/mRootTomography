function [rho, c] = rt_pinv(M, p, r)
%rt_PINV Summary of this function goes here
%   Detailed explanation goes here

s = size(M,1);
B = rt_meas_matrix(M);

rho = reshape(B\p, s, s);
tr = trace(rho);
[U,D] = svd(rho);
U = U(:,1:r);
D = D(1:r,1:r);
D = D / trace(D) * tr;
c = U*sqrt(D);
rho = c*c';

end

