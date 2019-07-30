function [rho, B] = rt_dm_pinv(k, M, n, r)
%rt_PINV Summary of this function goes here
%   Detailed explanation goes here

s = size(M,1);
B = rt_dm_B(M);

rho = pinv(B.*repmat(n,1,s^2))*k;
rho = reshape(rho, s, s);
[U,S] = svd(rho);
S((r+1):end, (r+1):end) = 0;
rho = U*S*U';
rho = rho/trace(rho);

end

