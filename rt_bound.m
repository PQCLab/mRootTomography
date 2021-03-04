function d = rt_bound(dm, proto, nshots, objType, varargin)
%rt_CHI_THEORY TODO
p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'dm');
addRequired(p, 'proto');
addRequired(p, 'nshots');
addRequired(p, 'objType');
addParameter(p, 'rank', 'dm');
addParameter(p, 'tracePreserving', true);
parse(p, dm, proto, nshots, objType, varargin{:});
opt = p.Results;

dim = size(dm, 1);
dm = dm / trace(dm);
H = rt_infomatrix(dm, proto, nshots, opt.objType, 'rank', opt.rank);

Constraints = [];
[Uh, Sh] = eig(H);
h = diag(Sh);

if strcmpi(opt.objType, 'state') || (strcmpi(opt.objType, 'process') && ~opt.tracePreserving)
    % Normalization constraint
    if ischar(opt.rank) && strcmp(opt.rank,'dm')
        opt.rank = rank(dm);
    end
    c = rt_purify(dm, opt.rank);
    Constraints = horzcat(Constraints, [real(c(:)); imag(c(:))]);
elseif strcmpi(opt.objType, 'process')
    % Trace preserving constraint
    if ischar(opt.rank) && strcmp(opt.rank, 'dm')
        opt.rank = rank(dm);
    end
    E = rt_process_reform(dm, 'chi2kraus', opt.rank);
    ds2r = dim * opt.rank;
    ConsTP = zeros(2*ds2r, dim);
    dim_sq = sqrt(dim);
    for m = 1:dim_sq
        for n = 1:dim_sq
            ind = (m-1)*dim_sq + n;

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
end

% Phase insensitivity constraint
tol = max(h)*1e-10;
ind0 = find(h < tol);
if length(ind0) > opt.rank^2
   warning('Information matrix has more than r^2 zero eigenvalues');
end
h(ind0) = tol;
ind0 = ind0(1:opt.rank^2);
Constraints = horzcat(Constraints, Uh(:,ind0));

% Finding d
L = null(Constraints')'*Uh*diag(1./sqrt(h));
d = svd(L).^2;

if strcmpi(opt.objType, 'process')
    d = d / sqrt(dim);
end

end

