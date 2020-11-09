function [chi, rinfo] = rt_chi_reconstruct(clicks, proto, varargin)
%rt_CHI_RECONSTRUCT Performs the reconstruction of chi-matrix, normalized
%to unity. The reconstruction is by monotone proximal gradient descend based on the paper
%https://papers.nips.cc/paper/5728-accelerated-proximal-gradient-methods-for-nonconvex-programming.pdf
% TODO

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'clicks');
addRequired(p, 'proto');
addOptional(p, 'nshots', 'sum');
addParameter(p, 'rank', 'auto');
addParameter(p, 'normalize', true);
addParameter(p, 'init', 'pinv');
addParameter(p, 'pinvOnly', false);
addParameter(p, 'tracePreserving', true);
addParameter(p, 'significanceLevel', 0.05);
addParameter(p, 'tol', 1e-8, @(x)x>0);
addParameter(p, 'maxIter', 1e6, @(x)x>0);
addParameter(p, 'display', false);
parse(p,clicks,proto,varargin{:});
op = p.Results;

% Recursive method for automatic rank estimation
d2 = size(proto{1},1);
d = sqrt(d2);
if ischar(op.rank) && strcmpi(op.rank, 'auto')
    [chi_r, rinfo_r] = deal(cell(1,d2));
    for i = 1:d2
        r = i;
        [chi_r{r}, rinfo_r{r}] = rt_chi_reconstruct(clicks, proto, varargin{:}, 'Rank', r);
        if isnan(rinfo_r{r}.pval) || rinfo_r{r}.pval >= op.significanceLevel
            break;
        elseif r > 1 && rinfo_r{r}.pval < rinfo_r{r-1}.pval
            r = r - 1;
            break;
        end
    end
    chi = chi_r{r};
    rinfo = rinfo_r{r};
    rinfo.chi_r = chi_r(1:r);
    rinfo.info_r = rinfo_r(1:r);
    return;
end

if ~op.tracePreserving
    [chi, rinfo] = rt_dm_reconstruct(clicks, proto, varargin{:}, 'Rank', op.rank);
    if op.normalize
        chi = chi * sqrt(d2);
    end
    return;
end

[proto,nshots] = rt_proto_check(proto, op.nshots, clicks);
if op.rank < 1 || op.rank > d2
    error('Process matrix rank should be between 1 and Hilbert space dimension');
end

ex = rt_experiment(d, 'poly');
ex = ex.set_data('proto', proto, 'nshots', nshots, 'clicks', clicks);

if strcmpi(op.init,'pinv') || op.pinvOnly
    p_est = ex.get_field('vec_clicks') ./ ex.get_field('vec_nshots');
    [~, e] = rt_pinv(cat(3, proto{:}), p_est, op.rank);
end

if ~strcmpi(op.init,'pinv')
    e = rt_purify(op.init, op.rank);
end

e = project_tp(e);

rinfo.iter = 0;
if ~op.pinvOnly
    optim = rt_optimizer('proximal_descend');
    optim = optim.set_options('display', op.display, 'tol', op.tol, 'max_iter', op.maxIter);
    [e, info] = optim.run(e, ...
        @(e) ex.get_logL_sq(e), ...                   %% log-likelihood
        @(e) ex.get_dlogL_sq(e), ...                  %% log-likelihood gradient
        @(e) project_tp(e / sqrt(trace(e'*e)/d)), ... %% proximal operation
        sum(nshots) ...                               %% Lipschitz constant
        );
    rinfo.iter = info.iter;
end

chi = e*e';
[rinfo.pval, rinfo.chi2, rinfo.df, rinfo.n_observed, rinfo.n_expected] = rt_significance(chi, clicks, proto, nshots, 'Rank', op.rank, 'isProcess', true);
rinfo.rank = op.rank;
if op.normalize
    chi = chi / trace(chi) * d;
end

end

function e = project_tp(e)
    [d2, r] = size(e);
    tr = trace(e*e');
    d = sqrt(d2);
    Ue = transpose(reshape(permute(reshape(e,[d,d,r]),[2,1,3]),d,[]));
    [U,~,V] = svd(Ue, 'econ'); Ue = U*V';
    e = reshape(permute(reshape(transpose(Ue),[d,d,r]),[2,1,3]),[d2,r]);
    e = e / sqrt(trace(e*e')/tr);
end
