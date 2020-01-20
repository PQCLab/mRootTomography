function d = rt_chi_theory(chi, proto, nshots, varargin)
%rt_CHI_THEORY TODO
p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'chi');
addRequired(p, 'proto');
addRequired(p, 'nshots');
addOptional(p, 'rank', 'dm');
addParameter(p, 'tracePreserving', true);
parse(p,chi,proto,nshots,varargin{:});
opt = p.Results;

s2 = size(chi,1);
s = sqrt(s2);
if ~opt.tracePreserving
    d = rt_dm_theory(chi, proto, nshots, varargin{:})/s;
    return;
end

chi = chi / trace(chi);
proto = rt_proto_check(proto, nshots);
H = rt_infomatrix(chi, proto, nshots, varargin{:});

Constraints = [];
[Uh,Sh] = eig(H);
h = diag(Sh);

% Trace preserving constraint
if ischar(opt.rank) && strcmp(opt.rank,'dm')
    opt.rank = rank(chi);
end
E = rt_root2kraus(rt_purify(chi,opt.rank));
ds2r = s2*opt.rank;
ConsTP = zeros(2*ds2r,s2);
for m = 1:s
    for n = 1:s
        ind = (m-1)*s + n;
        
        Deriv = zeros(size(E));
        Deriv(:,n,:) = Deriv(:,n,:) + E(:,m,:);
        Deriv(:,m,:) = Deriv(:,m,:) + conj(E(:,n,:));
        ConsTP(1:ds2r,ind) = Deriv(:);
        
        Deriv = zeros(size(E));
        Deriv(:,n,:) = Deriv(:,n,:) - 1j*E(:,m,:);
        Deriv(:,m,:) = Deriv(:,m,:) + 1j*conj(E(:,n,:));
        ConsTP((ds2r+1):end,ind) = Deriv(:);
    end
end
Constraints = horzcat(Constraints, ConsTP);

% Phase insensitivity constraint
tol = 0.01;
ind0 = find(h < tol);
if length(ind0) > opt.rank^2
   warning('Information matrix has more than r^2 zero eigenvalues');
end
h(ind0) = tol;
ind0 = ind0(1:opt.rank^2);
Constraints = horzcat(Constraints, Uh(:,ind0));

% Finding d
L = null(Constraints')'*Uh*diag(1./sqrt(2*h));
d = svd(L).^2 / s;

end

