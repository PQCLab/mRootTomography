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
addParameter(p, 'tracePreserving', true);
addParameter(p, 'significanceLevel', 0.05);
addParameter(p, 'display', 'none');
parse(p,clicks,proto,varargin{:});
opt = p.Results;

% Recursive method for automatic rank estimation
s = size(proto{1},1);
if ischar(opt.rank) && strcmpi(opt.rank, 'auto')
    [chi_r, rinfo_r] = deal(cell(1,s));
    for i = 1:s
        r = i;
        [chi_r{r}, rinfo_r{r}] = rt_chi_reconstruct(clicks, proto, varargin{:}, 'Rank', r);
        if isnan(rinfo_r{r}.pval) || rinfo_r{r}.pval >= opt.significanceLevel
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

if ~opt.tracePreserving
    [chi, rinfo] = rt_dm_reconstruct(clicks, proto, varargin{:}, 'Rank', opt.rank);
    if opt.normalize
        chi = chi * sqrt(s);
    end
    return;
end

nshots = opt.nshots;
[proto,nshots] = rt_proto_check(proto,nshots);

% Params
dispfreq = 50*1;
eps = 1e-10;
N_max = 100000;

if dispfreq > 0
    h = rt_fprintreplace('Initialization...');
end

% First approximation by pseudo-inverse
[M, n, k] = rt_data_join(proto, nshots, clicks);
[~, e] = rt_pinv(M, k./n, opt.rank);
e = project_tp(e);
B = rt_meas_matrix(M);

% Accelerated projective gradient descend
L = sum(n)*10;
x0 = e; x1 = e; z1 = e;
t0 = 0; t1 = 1;
alp = 1/L;
F_x1 = 0;
for i = 1:N_max
    y1 = x1 + t0/t1*(z1-x1) + (t0-1)/t1*(x1-x0);
    z2 = project_tp(y1-alp*get_grad(y1, k, B, n, s));
    v2 = project_tp(x1-alp*get_grad(x1, k, B, n, s));
    F_z2 = get_func(z2,k,n,B);
    F_v2 = get_func(v2,k,n,B);
    if F_z2 <= F_v2
        x2 = z2;
        F_x2 = F_z2;
    else
        x2 = v2;
        F_x2 = F_v2;
    end
    
    x0 = x1;
    x1 = x2;
    z1 = z2;
    
    t0 = t1;
    t1 = (sqrt(4*t0^2+1)+1)/2;
    
    dx = norm(x1-x0);
    dF = abs(F_x2-F_x1);
    F_x1 = F_x2;
    stopIter = (dx < eps);
    
    if dispfreq > 0 && (mod(i,dispfreq) == 0 || i == 1 || stopIter)
        h = rt_fprintreplace(sprintf('r = %d | dx = %.4e | dF = %.4e | F = %.4e | %d', opt.rank, dx, dF, F_x1, i), h);
    end
    if stopIter
        break;
    end
end
if dispfreq > 0
    fprintf('\n');
end

chi = x1*x1';
[rinfo.pval, rinfo.chi2, rinfo.df, rinfo.n_observed, rinfo.n_expected] = rt_significance(chi, clicks, proto, nshots, 'Rank', opt.rank, 'isProcess', true);
rinfo.iter = i;
rinfo.rank = opt.rank;
if opt.normalize
    chi = chi / trace(chi) * sqrt(s);
end

end

function F = get_func(e, k, n, B)
pn = get_p(e,B).*n;
F = -sum(k.*log(pn)-pn);
end

function D = get_grad(e, k, B, n, s)
a = n - k./get_p(e,B);
J = reshape(B'*a, s, s);
D = 2*J*e;
end

function p = get_p(e, B)
p = real(B * reshape(e*e',[],1));
end

function e = project_tp(e)
[s2, r] = size(e);
tr = trace(e*e');
s = sqrt(s2);
Ue = transpose(reshape(permute(reshape(e,[s,s,r]),[2,1,3]),s,[]));
[U,~,V] = svd(Ue, 'econ'); Ue = U*V';
e = reshape(permute(reshape(transpose(Ue),[s,s,r]),[2,1,3]),[s2,r]);
e = e / sqrt(trace(e*e')/tr);
end
