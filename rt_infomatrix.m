function h = rt_infomatrix(dm, ex, varargin)
% RT_INFOMATRIX Calculates the complete Fisher information matrix by the
% density matrix and the measurements protocol
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
opt.rank = 'dm';
for ja = 1:2:length(varargin)
    opt.(lower(varargin{ja})) = varargin{ja + 1};
end

dm_norm = trace(dm);
if ischar(opt.rank) && strcmp(opt.rank,'dm')
    opt.rank = rank(dm);
end
psi = rt_purify(dm, opt.rank);

% Find close state with no zeros probabilities
pTol = 1e-8;
Ntries = 100;
for i = 1:Ntries
    prob = ex.get_probs_sq(psi);
    if any(prob < pTol | prob > 1 - pTol)
        if i == Ntries
            warning('RT:NonSingulatState', 'Failed to find non-singular state');
        else
            psi = psi + (randn(size(psi)) + 1j*randn(size(psi))) * sqrt(pTol);
            psi = psi / sqrt(trace(psi'*psi) / dm_norm);
        end
    else
        break;
    end
end

% Calculate fisher information
h = 0;
operators = cat(3, ex.proto{:});
nshots = ex.get_field('vec_nshots');
for j = 1:size(operators, 3)
    a = reshape(operators(:,:,j) * psi, [], 1);
    a = [real(a); imag(a)]; 
    h = h + ex.stat().fisher_information(nshots(j), prob(j)) * (a * a');
end
h = 4 * h;

end

