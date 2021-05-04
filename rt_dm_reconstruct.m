function [dm, rinfo] = rt_dm_reconstruct(ex, varargin)
% RT_DM_RECONSTRUCT Reconstruct the quantum state density matrix by the
% results of a set of complementary measurements
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'ex', @(ex) isa(ex, 'rt_experiment'));
addParameter(p, 'rank', 'auto');
addParameter(p, 'significanceLevel', 0.05, @(x)x>0&&x<1);
addParameter(p, 'getStats', false);
addParameter(p, 'init', 'pinv');
addParameter(p, 'regCoeff', 0.5, @(x)x>0&&x<=1);
addParameter(p, 'tol', 1e-8, @(x)x>0);
addParameter(p, 'maxIter', 1e6, @(x)x>0);
addParameter(p, 'display', false);
parse(p, ex, varargin{:});
op = p.Results;

dim = ex.dim;
if ischar(op.rank) && strcmpi(op.rank, 'auto')
    optim = rt_optimizer('auto_rank');
    optim.set_options('display', op.display, 'sl', op.significanceLevel);
    [data, data_r] = optim.run(dim, @(r) rank_fun(ex, varargin, r));
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

if strcmpi(op.init, 'pinv')
    p_est = ex.get_field('vec_clicks') ./ ex.get_field('vec_nshots');
    [~, psi] = rt_pinv(ex.get_field('vec_proto'), p_est, op.rank);
else
    psi = rt_purify(op.init, op.rank);
end

optim = rt_optimizer('fixed_point');
optim = optim.set_options('display', op.display, 'tol', op.tol, 'max_iter', op.maxIter, 'reg_coeff', op.regCoeff);
mu_inv = 1 / ex.logL_eq_mu();
foptim = @(psi) mu_inv * ex.logL_eq_jmat_dm(psi * psi' / trace(psi' * psi)) * psi;
[psi, optim_info] = optim.run(psi, foptim);
rinfo.optimizer = optim;
rinfo.iter = optim_info.iter;

dm = psi*psi';
dm = dm / trace(dm);

rinfo.rank = op.rank;
if op.getStats
    rinfo.chi2 = ex.chi2_dm(dm);
    rinfo.df = ex.deg_f_rank(op.rank);
    rinfo.pval = rt_pval(rinfo.chi2, rinfo.df);
end

end

function data = rank_fun(ex, args, r)
    [dm, data] = rt_dm_reconstruct(ex, args{:}, 'GetStats', true, 'Rank', r);
    data.dm = dm;
end