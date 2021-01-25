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
addParameter(p, 'init', 'pinv');
addParameter(p, 'pinvOnly', false);
addParameter(p, 'getStats', false);
addParameter(p, 'tracePreserving', true);
addParameter(p, 'significanceLevel', 0.05);
addParameter(p, 'tol', 1e-8, @(x)x>0);
addParameter(p, 'maxIter', 1e6, @(x)x>0);
addParameter(p, 'display', false);
parse(p,clicks,proto,varargin{:});
op = p.Results;

[proto, nshots] = rt_proto_check(proto, op.nshots, clicks);
dim2 = size(proto{1}, 1);
dim = sqrt(dim2);

if ischar(op.rank) && strcmpi(op.rank, 'auto')
    [data, ~, data_r] = rt_optimize_rank(dim2, op.significanceLevel, op.display, @(r) rank_fun(clicks, proto, varargin, r));
    chi = data.chi;
    rinfo = rmfield(data, 'chi');
    rinfo.data_r = data_r;
    return;
elseif ischar(op.rank) && strcmpi(op.rank, 'full')
    op.rank = dim2;
end

if op.rank < 1 || op.rank > dim2
    error('RT:RankValue', 'Density matrix rank should be between 1 and squared Hilbert space dimension');
end

if ~op.tracePreserving
    [chi, rinfo] = rt_dm_reconstruct(clicks, proto, varargin{:}, 'Rank', op.rank);
    chi = chi * sqrt(dim);
    return;
end

ex = rt_experiment(dim2, 'poly', 'process');
ex.set_data('proto', proto, 'nshots', nshots, 'clicks', clicks);

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
    optim.set_options('display', op.display, 'tol', op.tol, 'max_iter', op.maxIter);
    [e, info] = optim.run(e, ...
        @(e) ex.get_logL_sq(e), ...                     %% log-likelihood
        @(e) ex.get_dlogL_sq(e), ...                    %% log-likelihood gradient
        @(e) project_tp(e / sqrt(trace(e'*e)/dim)), ... %% proximal operation
        sum(nshots) ...                                 %% Lipschitz constant
        );
    rinfo.optimizer = optim;
    rinfo.iter = info.iter;
end
chi = e*e';
chi = chi / trace(chi) * dim;

rinfo.rank = op.rank;
rinfo.experiment = ex;
if op.getStats
    rinfo.chi2 = ex.get_chi2_dm(chi);
    rinfo.df = ex.get_df(op.rank);
    rinfo.pval = rt_pval(rinfo.chi2, rinfo.df);
end

end

function e = project_tp(e)
    [d2, r] = size(e);
    tr = trace(e*e');
    d = sqrt(d2);
    ue = transpose(reshape(permute(reshape(e,[d,d,r]),[2,1,3]),d,[]));
    [u,~,v] = svd(ue, 'econ');
    ue = u*v';
    e = reshape(permute(reshape(transpose(ue),[d,d,r]),[2,1,3]),[d2,r]);
    e = e / sqrt(trace(e*e')/tr);
end

function data = rank_fun(clicks, proto, args, r)
    [chi, data] = rt_chi_reconstruct(clicks, proto, args{:}, 'GetStats', true, 'Rank', r);
    data.chi = chi;
end
