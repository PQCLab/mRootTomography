function [rho, rinfo] = rt_dm_reconstruct(clicks, proto, varargin)
%RT_DM_RECONSTRUCT TODO
%   TODO

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'clicks');
addRequired(p, 'proto');
addOptional(p, 'nshots', 'sum');
addParameter(p, 'rank', 'auto');
addParameter(p, 'normalize', true);
addParameter(p, 'pinvOnly', false);
addParameter(p, 'significanceLevel', 0.5);
addParameter(p, 'display', 'none');
parse(p,clicks,proto,varargin{:});
opt = p.Results;

% Recursive method for automatic rank estimation
s = size(proto{1},1);
if ischar(opt.rank) && strcmpi(opt.rank, 'auto')
    [rho_r, rinfo_r] = deal(cell(1,s));
    for i = 1:s
        r = i;
        [rho_r{r}, rinfo_r{r}] = rt_dm_reconstruct(clicks, proto, opt.nshots, 'Rank', r);
        if isnan(rinfo_r{r}.pval) || rinfo_r{r}.pval >= opt.significanceLevel
            break;
        elseif r > 1 && rinfo_r{r}.pval < rinfo_r{r-1}.pval
            r = r - 1;
            break;
        end
    end
    rho = rho_r{r};
    rinfo = rinfo_r{r};
    rinfo.rho_r = rho_r(1:r);
    rinfo.info_r = rinfo_r(1:r);
    return;
end

if ~iscell(proto)
    proto = num2cell(proto,[1 2]);
end
rt_proto_check(proto, opt.nshots);

% Params
eps = 1e-10;
N_max = 100000;
alp = 0.5;

% First approximation by pseudo-inverse
[M, n, k] = rt_data_join(proto, opt.nshots, clicks);
[~, c] = rt_pinv(M, k./n, opt.rank);
B = rt_meas_matrix(M);

i = 0;
if ~opt.pinvOnly
    s = size(c,1);
    Ir = inv(reshape(B'*n, s, s));

    % Iterations
    for i = 1:N_max
        cp = c;
        c = alp*get_new_value(c, k, B, Ir, s) + (1-alp)*cp;
        dc = norm(cp-c);
        if (dc < eps)
            break;
        end
    end
end

rho = c*c';
[rinfo.pval, rinfo.chi2, rinfo.df] = rt_significance(rho, proto, opt.nshots, clicks, opt.rank, false);
rinfo.iter = i;
rinfo.rank = opt.rank;
if opt.normalize
    rho = rho / trace(rho);
end

end

function c = get_new_value(c, k, B, Ir, s)
rho = c*c';
a = k./(B*rho(:));
J = reshape(B'*a, s, s);
c = Ir*J*c;
end