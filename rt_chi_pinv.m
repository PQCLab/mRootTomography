function [chi, B, kv] = rt_chi_pinv(k, M, n, r)
%rt_PINV Reconstruction by pseudo-inverse
%   INPUT: same as for rt_chi_reconstruct
%
%   OUTPUT:
%   chi - chi-matrix (normalized to unity)
%   B - measurement matrix
%   kv - vectorized measurement counts

s = size(M{1},1);

B = rt_chi_B(M,n);
kv = k(:);

chi = pinv(B)*kv;
chi = reshape(chi, s, s);
[U,S] = svd(chi);
S((r+1):end, (r+1):end) = 0;
chi = U*S*U';
chi = chi/trace(chi);

end

