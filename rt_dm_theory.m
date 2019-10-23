function d = rt_dm_theory(dm, proto, nshots, r)
%RT_DM_THEORY Summary of this function goes here TODO
%   Detailed explanation goes here TODO

dm = dm / trace(dm);
if nargin < 4 || strcmp(r,'auto')
    r = rank(dm);
end
proto = rt_proto_check(proto, nshots);
H = rt_infomatrix(dm, proto, nshots, r);

% Normalization constraint
c = rt_purify(dm,r);
[M, n] = rt_data_join(proto, nshots);
R = 2*sum(bsxfun(@times, M, reshape(n,1,1,[])), 3) * c;
R = [real(R(:)); imag(R(:))];
U = null(R');
H = U'*H*U;

% Finding d
h = svd(H);
ind0 = find(h < 1e-4);
if length(ind0) > r^2
   warning('Information matrix has more than r^2 zero eigenvalues. Extra zeros are ignored');
end
h(ind0) = [];
d = 1./(2*h);

end

