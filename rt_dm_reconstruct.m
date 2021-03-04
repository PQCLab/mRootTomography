function [dm, rinfo] = rt_dm_reconstruct(dim, clicks, proto, varargin)
%RT_DM_RECONSTRUCT Performs the quantum state reconstruction using root
%approach and maximum likelihood. Detail explanation of the algorithm could
%be found in the paper [Borganov Yu. I., JETP 108(6) 2009] 
%
%   dm = rt_dm_reconstruct(clicks,proto) reconstructs density matrix by
%   the set of clicks obtained in each measurement experiment defined by
%   protocol proto. clicks{j} - column-vector with the numbers of observed
%   events in j-th experiment. proto{j} - 3D-array of measurement operators.
%   size(proto{j},3) should be equal to length(click{j}) for every j.
%
%   dm = rt_dm_reconstruct(clicks,proto,nshots) specifies the row-vector
%   nshots - row-vector of the numbers of measurements. nshots(j) - number
%   of measurements for j-th experiment. length(nshots) should be equal to
%   length(proto) and length(clicks). If nshots = 'sum' then for j-th
%   experiment number of measurements is sum(clicks{j}). Default: 'sum'
%
%   dm = rt_dm_reconstruct( ___ ,Name,Value) specifies additional parameters
%   for reconstruction. For example, rt_dm_reconstruct(clicks,proto,'Rank',1)
%   reconstructs a pure (rank-1) quantum state. Available parameters:
%       • Rank (integer or sting) - rank of the quantum state model. 'auto'
%       indicates that rank should be calculated automatically using
%       chi-squared test, 'full' indicates full-rank reonstruction. Default: 'auto'
%       • Normalize (boolean) - normalize output density matrix. Default: true
%       • Init (string or dm) - initial guess of the density matrix. Value
%       'pinv' stands for the initial guess by pseudo-inversion.
%       Default: 'pinv'
%       • PinvOnly (boolean) - reconstruct density matrix by the
%       pseudo-inversion only. Default: false
%       • SignificanceLevel (float) - significance level for the
%       automatic rank definition (used when Rank equal to 'auto').
%       Default: 0.005
%       • Alpha (float) - convergence speed parameter. At each step old state
%       matrix C becomes Alpha*CN + (1-Alpha)*C, where CN is a new state matrix.
%       Default: 0.5
%       • Tol (float) - termination tolerance on state matrix. Default: 1e-8
%       • MaxIter (integer) - maximum number of iterations. Default: 1e6
%       • Display (bool) - display iterations. Default: false
%
%   [dm,rinfo] = rt_dm_reconstruct( ___ ) also returns additional
%   reconstruction informations
%
%OUTPUT:
%   dm - density matrix
%   rinfo - structure array with the following fields:
%       • iter - number of iteration taken to find likelihood maximum
%       • rank - rank of the quantum state model
%       • pval, chi2, df - p-value, chi-squared value and number of degrees
%        of freedom of the statistical chi-squared test
%       • dm_r, info_r - cell arrays of density matrices and reconstruction
%       information for every value of rank below optimal one (defined if
%       parameter 'Rank' is set to 'auto')
%
%Author: PQCLab
%Website: https://github.com/PQCLab/RootTomography
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