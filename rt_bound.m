function d = rt_bound(dm, ex, varargin)
% RT_BOUND Calculates the lower bound for the variances of quantum state or
% quantum process parameters estimator
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
opt.rank = 'dm';
for ja = 1:2:length(varargin)
    opt.(lower(varargin{ja})) = varargin{ja + 1};
end

dim = size(dm, 1);
H = rt_infomatrix(dm, ex, 'rank', opt.rank);

Constraints = [];
[Uh, Sh] = eig(H);
h = diag(Sh);

if ischar(opt.rank) && strcmp(opt.rank,'dm')
    opt.rank = rank(dm);
end

if strcmp(ex.obj_type, 'state')
    % Normalization constraint
    psi = rt_purify(dm, opt.rank);
    Constraints = horzcat(Constraints, [real(psi(:)); imag(psi(:))]);
elseif strcmpi(ex.obj_type, 'process')
    % Trace preserving constraint
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
d = d / trace(dm);

end

