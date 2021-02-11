function h = rt_infomatrix(dm, proto, nshots, varargin)
% RT_INFOMATRIX TODO
p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'dm');
addRequired(p, 'proto');
addRequired(p, 'nshots');
addOptional(p, 'rank', 'dm');
parse(p,dm,proto,nshots,varargin{:});
opt = p.Results;

if ischar(opt.rank) && strcmp(opt.rank,'dm')
    opt.rank = rank(dm);
end

dim = size(dm, 1);
c = rt_purify(dm, opt.rank);
ex = rt_experiment(dim, 'auto', 'state');
ex.set_data('proto', proto, 'nshots', nshots);

% Find close state with no zeros probabilities
pTol = 1e-10;
Ntries = 100;
for i = 1:Ntries
    prob = ex.get_probs_sq(c);
    if any(prob < pTol)
        if i == Ntries
            warning('RT:NonSingulatState', 'Failed to find non-singular state');
        else
            c = c + (randn(size(c)) + 1j*randn(size(c))) * sqrt(pTol);
            c = c / sqrt(trace(c'*c));
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
    a = reshape(operators(:,:,j) * c, [], 1);
    a = [real(a); imag(a)]; 
    h = h + nshots(j) * (a * a') / prob(j);
end
h = 2*h;

end

