function d = rt_dm_theory(dm, proto, nshots, varargin)
%RT_DM_THEORY TODO
p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'dm');
addRequired(p, 'proto');
addRequired(p, 'nshots');
addOptional(p, 'rank', 'dm');
parse(p,dm,proto,nshots,varargin{:});
opt = p.Results;

dm = dm / trace(dm);
proto = rt_proto_check(proto, nshots);
H = rt_infomatrix(dm, proto, nshots, varargin{:});

% Normalization constraint
if ischar(opt.rank) && strcmp(opt.rank,'dm')
    opt.rank = rank(dm);
end
c = rt_purify(dm,opt.rank);
[M, n] = rt_data_join(proto, nshots);
R = 2*sum(bsxfun(@times, M, reshape(n,1,1,[])), 3) * c;
R = [real(R(:)); imag(R(:))];
U = null(R');
H = U'*H*U;

% Finding d
h = svd(H);
ind0 = find(h < 1e-4);
if length(ind0) > opt.rank^2
   warning('Information matrix has more than r^2 zero eigenvalues. Extra zeros are ignored');
end
h(ind0) = [];
d = 1./(2*h);

end

