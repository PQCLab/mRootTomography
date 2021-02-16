function [chi, rinfo] = rt_chi_reconstruct(dim, clicks, proto, varargin)
%rt_CHI_RECONSTRUCT Performs the reconstruction of chi-matrix, normalized
%to unity. The reconstruction is by monotone proximal gradient descend based on the paper
%https://papers.nips.cc/paper/5728-accelerated-proximal-gradient-methods-for-nonconvex-programming.pdf
% TODO

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'dim');
addRequired(p, 'clicks');
addRequired(p, 'proto');
addOptional(p, 'nshots', 'sum');
addParameter(p, 'rank', 'auto');
addParameter(p, 'statType', 'auto');
addParameter(p, 'init', 'pinv');
addParameter(p, 'pinvOnly', false);
addParameter(p, 'getStats', false);
addParameter(p, 'tracePreserving', true);
addParameter(p, 'significanceLevel', 0.05);
addParameter(p, 'tol', 1e-8, @(x)x>0);
addParameter(p, 'maxIter', 1e6, @(x)x>0);
addParameter(p, 'display', false);
parse(p, dim, clicks, proto, varargin{:});
op = p.Results;

if ischar(op.rank) && strcmpi(op.rank, 'auto')
    [data, ~, data_r] = rt_optimize_rank(dim^2, op.significanceLevel, op.display, @(r) rank_fun(dim, clicks, proto, varargin, r));
    chi = data.chi;
    rinfo = rmfield(data, 'chi');
    rinfo.data_r = data_r;
    return;
elseif ischar(op.rank) && strcmpi(op.rank, 'full')
    op.rank = dim^2;
end

if op.rank < 1 || op.rank > dim^2
    error('RT:RankValue', 'Process matrix rank should be between 1 and squared Hilbert space dimension');
end

ex = rt_experiment(dim, 'process', op.statType);
ex.set_data('proto', op.proto, 'clicks', op.clicks);
if ischar(op.nshots) && strcmpi(op.nshots, 'sum')
    op.nshots = cellfun(@(kj) sum(kj), ex.clicks);
end
ex.set_data('nshots', op.nshots);

if op.tracePreserving
    if strcmpi(op.init, 'pinv') || op.pinvOnly
        p_est = ex.get_field('vec_clicks') ./ ex.get_field('vec_nshots');
        [~, e] = rt_pinv(cat(3, ex.proto{:}), p_est, op.rank);
    end
    if ~strcmpi(op.init, 'pinv')
        e = rt_purify(op.init, op.rank);
    end
    e = project_tp(e);
    rinfo.iter = 0;
    if ~op.pinvOnly
        optim = rt_optimizer('proximal_descend');
        optim.set_options('display', op.display, 'tol', op.tol, 'max_iter', op.maxIter);
        [e, info] = optim.run(e, ...
            @(e) ex.get_logL_sq(e), ...                         %% log-likelihood
            @(e) ex.get_dlogL_sq(e), ...                        %% log-likelihood gradient
            @(e) project_tp(e / sqrt(trace(e'*e) / dim)), ...   %% proximal operation
            10 * sum(op.nshots) ...                             %% Lipschitz constant
            );
        rinfo.optimizer = optim;
        rinfo.iter = info.iter;
    end
    chi = e*e';
    chi = chi / trace(chi) * dim;
else
    if strcmp(ex.stat_type, 'poly')
        warning('RT:PolyStatNonTPProcess', 'Polynomial statistics is incompatible with non trace preserving processes');
    end
    [chi, rinfo] = rt_dm_reconstruct(dim^2, ex.clicks, ex.proto, varargin{:}, 'statType', ex.stat_type, 'Rank', op.rank, 'getStats', false);
    chi = chi * dim;
    ex = rinfo.experiment;
end

rinfo.rank = op.rank;
rinfo.experiment = ex;
if op.getStats
    rinfo.chi2 = ex.get_chi2_dm(chi);
    rinfo.df = ex.get_df(op.rank);
    rinfo.pval = rt_pval(rinfo.chi2, rinfo.df);
end

end

function e = project_tp(e)
    [dim2, r] = size(e);
    tr = trace(e*e');
    dim = sqrt(dim2);
    ue = transpose(reshape(permute(reshape(e, [dim, dim, r]), [2, 1, 3]), dim, []));
    [u, ~, v] = svd(ue, 'econ');
    ue = u * v';
    e = reshape(permute(reshape(transpose(ue), [dim, dim, r]), [2, 1, 3]), [dim2, r]);
    e = e / sqrt(trace(e*e') / tr);
end

function data = rank_fun(dim, clicks, proto, args, r)
    [chi, data] = rt_chi_reconstruct(dim, clicks, proto, args{:}, 'GetStats', true, 'Rank', r);
    data.chi = chi;
end
