function [dm, rinfo] = rt_dm_reconstruct(dim, clicks, proto, varargin)
% RT_DM_RECONSTRUCT Reconstruct the quantum state density matrix by the
% results of a set of complementary measurements
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'dim');
addRequired(p, 'clicks');
addRequired(p, 'proto');
addOptional(p, 'nshots', 'sum');
addParameter(p, 'statType', 'auto');
addParameter(p, 'rank', 'auto');
addParameter(p, 'significanceLevel', 0.05, @(x)x>0&&x<1);
addParameter(p, 'getStats', false);
addParameter(p, 'init', 'pinv');
addParameter(p, 'regCoeff', 0.5, @(x)x>0&&x<=1);
addParameter(p, 'tol', 1e-8, @(x)x>0);
addParameter(p, 'maxIter', 1e6, @(x)x>0);
addParameter(p, 'display', false);
parse(p, dim, clicks, proto, varargin{:});
op = p.Results;

if ischar(op.rank) && strcmpi(op.rank, 'auto')
    optim = rt_optimizer('auto_rank');
    optim.set_options('display', op.display, 'sl', op.significanceLevel);
    [data, data_r] = optim.run(dim, @(r) rank_fun(dim, clicks, proto, varargin, r));
    dm = data.dm;
    rinfo = rmfield(data, 'dm');
    rinfo.data_r = data_r;
    return;
elseif ischar(op.rank) && strcmpi(op.rank, 'full')
    op.rank = dim;
end

if op.rank < 1 || op.rank > dim
    error('RT:RankValue', 'Density matrix rank should be between 1 and Hilbert space dimension');
end

ex = rt_experiment(dim, 'state', op.statType);
ex.set_data('proto', op.proto, 'clicks', op.clicks);
if strcmpi(op.nshots, 'sum')
    op.nshots = cellfun(@(kj) sum(kj), ex.clicks);
end
ex.set_data('nshots', op.nshots);

if strcmpi(op.init, 'pinv')
    p_est = ex.get_field('vec_clicks') ./ ex.get_field('vec_nshots');
    [~, c] = rt_pinv(cat(3, ex.proto{:}), p_est, op.rank);
else
    c = rt_purify(op.init, op.rank);
end

optim = rt_optimizer('fixed_point');
optim = optim.set_options('display', op.display, 'tol', op.tol, 'max_iter', op.maxIter, 'reg_coeff', op.regCoeff);
Ir = inv(reshape(2 * (ex.get_field('vec_nshots')' * ex.get_field('vec_proto'))', size(c, 1), []));
if strcmp(ex.stat_type, 'poly')
    foptim = @(c) Ir * ex.get_dlogL_sq(c);
elseif strcmp(ex.stat_type, 'poiss')
    foptim = @(c) Ir * ex.get_dlogL_sq(c) + c;
else
    error('RT:StatsTypeRec', 'Only `poly`, `poiss` and `bino` statistics are currently supported in rt_dm_reconstruct');
end
[c, optim_info] = optim.run(c, foptim);
rinfo.optimizer = optim;
rinfo.iter = optim_info.iter;

dm = c*c';
dm = dm / trace(dm);

rinfo.rank = op.rank;
rinfo.experiment = ex;
if op.getStats
    rinfo.chi2 = ex.get_chi2_dm(dm);
    rinfo.df = ex.get_df(op.rank);
    rinfo.pval = rt_pval(rinfo.chi2, rinfo.df);
end

end

function data = rank_fun(dim, clicks, proto, args, r)
    [dm, data] = rt_dm_reconstruct(dim, clicks, proto, args{:}, 'GetStats', true, 'Rank', r);
    data.dm = dm;
end