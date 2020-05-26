function [dm, rinfo] = rt_dm_reconstruct(clicks, proto, varargin)
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
addRequired(p, 'clicks');
addRequired(p, 'proto');
addOptional(p, 'nshots', 'sum');
addParameter(p, 'rank', 'auto');
addParameter(p, 'normalize', true);
addParameter(p, 'init', 'pinv');
addParameter(p, 'pinvOnly', false);
addParameter(p, 'significanceLevel', 0.05, @(x)x>0&&x<1);
addParameter(p, 'alpha', 0.5, @(x)x>0&&x<=1);
addParameter(p, 'tol', 1e-8, @(x)x>0);
addParameter(p, 'maxIter', 1e6, @(x)x>0);
addParameter(p, 'display', false);
parse(p,clicks,proto,varargin{:});
opt = p.Results;

% Recursive method for automatic rank estimation
s = size(proto{1},1);
if ischar(opt.rank) && strcmpi(opt.rank, 'auto')
    if opt.display
        fprintf('=== Automatic rank estimation ===\n');
    end
    pvalRed = false;
    [dm_r, rinfo_r] = deal(cell(1,s));
    for i = 1:s
        r = i;
        if opt.display
            fprintf('Try rank %d\n', r);
        end
        [dm_r{r}, rinfo_r{r}] = rt_dm_reconstruct(clicks, proto, varargin{:}, 'Rank', r);
        if isnan(rinfo_r{r}.pval) || rinfo_r{r}.pval >= opt.significanceLevel
            break;
        elseif r > 1 && rinfo_r{r}.pval < rinfo_r{r-1}.pval
            pvalRed = true;
            r = r - 1;
            break;
        end
    end
    if opt.display
        if rinfo_r{r}.pval >= opt.significanceLevel
            fprintf('Rank %d is statistically significant at significance level %s. Procedure terminated.\n', r, num2str(opt.significanceLevel));
        elseif pvalRed
            fprintf('P-value is maximal (%s) for rank %d. Procedure terminated.\n', num2str(rinfo_r{r}.pval), r);
        else
            fprintf('Failed to determine optimal rank. Maximal rank %d is taken.\n', r);
        end
    end
    dm = dm_r{r};
    rinfo = rinfo_r{r};
    rinfo.dm_r = dm_r(1:r);
    rinfo.info_r = rinfo_r(1:r);
    return;
elseif ischar(opt.rank) && strcmpi(opt.rank, 'full')
    opt.rank = s;
end

[proto,nshots] = rt_proto_check(proto,opt.nshots,clicks);
if opt.rank < 1 || opt.rank > size(proto{1},1)
    error('Density matrix rank should be between 1 and Hilbert space dimension');
end

% Params
dispfreq = 50;
eps = opt.tol;
N_max = opt.maxIter;
alp = opt.alpha;

[M, n, k] = rt_data_join(proto, nshots, clicks);
if strcmpi(opt.init,'pinv') || opt.pinvOnly
    if opt.display
        h = rt_fprintreplace('Pseudo-inversion...');
    end
    [~, c] = rt_pinv(M, k./n, opt.rank);
end

if ~strcmpi(opt.init,'pinv')
    c = rt_purify(opt.init, opt.rank);
end

i = 0;
if ~opt.pinvOnly
    s = size(c,1);
    B = rt_meas_matrix(M);
    Ir = inv(reshape(B'*n, s, s));
    for i = 1:N_max
        cp = c;
        c = (1-alp)*get_new_value(c, k, B, Ir, s) + alp*cp;
        dc = norm(cp-c);
        stopIter = (dc < eps);
        
        if opt.display && (mod(i,dispfreq) == 0 || i == 1 || stopIter)
            h = rt_fprintreplace(sprintf('Iteration %d \t\t Difference %.4e', i, dc), h);
        end
        if stopIter
            break;
        end
    end
end
if opt.display
    fprintf('\n');
end

dm = c*c';
rinfo.iter = i;
rinfo.rank = opt.rank;
[rinfo.pval, rinfo.chi2, rinfo.df, rinfo.n_observed, rinfo.n_expected] = rt_significance(dm, clicks, proto, nshots, 'Rank', opt.rank);
if opt.normalize
    dm = dm / trace(dm);
end

end

function c = get_new_value(c, k, B, Ir, s)
rho = c*c';
p = real(B*rho(:));
p(p<1e-15) = 1e-15;
a = k./p;
J = reshape(B'*a, s, s);
c = Ir*J*c;
end