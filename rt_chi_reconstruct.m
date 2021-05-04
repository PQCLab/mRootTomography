function [chi, rinfo] = rt_chi_reconstruct(ex, varargin)
% RT_CHI_RECONSTRUCT Reconstruct the quantum process matrix by the results
% of a set of complementary measurements
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'ex', @(ex) isa(ex, 'rt_experiment'));
addParameter(p, 'statType', 'auto');
addParameter(p, 'rank', 'auto');
addParameter(p, 'significanceLevel', 0.05);
addParameter(p, 'getStats', false);
addParameter(p, 'init', 'pinv');
addParameter(p, 'LipsConst', 'ntot');
addParameter(p, 'tol', 1e-8, @(x)x>0);
addParameter(p, 'maxIter', 1e6, @(x)x>0);
addParameter(p, 'display', false);
parse(p, ex, varargin{:});
op = p.Results;

dim = ex.dim;
if ischar(op.rank) && strcmpi(op.rank, 'auto')
    optim = rt_optimizer('auto_rank');
    optim.set_options('display', op.display, 'sl', op.significanceLevel);
    [data, data_r] = optim.run(dim^2, @(r) rank_fun(ex, varargin, r));
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

if strcmpi(op.init, 'pinv')
    p_est = ex.get_field('vec_clicks') ./ ex.get_field('vec_nshots');
    [~, e] = rt_pinv(ex.get_field('vec_proto'), p_est, op.rank);
else
    e = rt_purify(op.init, op.rank);
end
e = e / sqrt(trace(e'*e) / dim);
e = project_tp(e);

if strcmpi(op.LipsConst, 'ntot')
    op.LipsConst = sum(ex.get_field('nshots'));
end

optim = rt_optimizer('proximal_ascend');
optim.set_options('display', op.display, 'tol', op.tol, 'max_iter', op.maxIter);
[e, info] = optim.run(e, ...
    @(e) ex.logL_sq(e), ...                             %% log-likelihood
    @(e) ex.dlogL_sq(e), ...                            %% log-likelihood gradient
    @(e) project_tp(e / sqrt(trace(e'*e) / dim)), ...   %% proximal operation
    op.LipsConst ...                                    %% Lipschitz constant
    );
rinfo.optimizer = optim;
rinfo.iter = info.iter;
chi = e*e';
chi = chi / trace(chi) * dim;

rinfo.rank = op.rank;
if op.getStats
    rinfo.chi2 = ex.chi2_dm(chi);
    rinfo.df = ex.deg_f_rank(op.rank);
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

function data = rank_fun(ex, args, r)
    [chi, data] = rt_chi_reconstruct(ex, args{:}, 'GetStats', true, 'Rank', r);
    data.chi = chi;
end