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

Constraints = [];
[Uh,Sh] = eig(H);
h = diag(Sh);

% Normalization constraint
if ischar(opt.rank) && strcmp(opt.rank,'dm')
    opt.rank = rank(dm);
end
c = rt_purify(dm,opt.rank);
Constraints = horzcat(Constraints, [real(c(:)); imag(c(:))]);

% Phase insensitivity constraint
tol = max(h)*1e-10;
ind0 = find(h < tol);
if length(ind0) > opt.rank^2
   warning('Information matrix has more than r^2 zero eigenvalues');
end
h(ind0) = tol;
ind0 = ind0(1:opt.rank^2);
Constraints = horzcat(Constraints, Uh(:,ind0));

% Variances
L = null(Constraints')'*Uh*diag(1./sqrt(2*h));
d = svd(L).^2;

end

