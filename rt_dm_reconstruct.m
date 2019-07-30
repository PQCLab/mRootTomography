function [rho, rinfo] = rt_dm_reconstruct(k0, M0, n0, r)
%RECONSTRUCT Summary of this function goes here
%   Detailed explanation goes here

% Recursive method for automatic rank estimation
s = size(k0,1);
if (ischar(r) && strcmpi(r, 'auto')) || r == 0
    [rho, rinfo] = rt_dm_reconstruct(k0,M0,n0,1);
    for r = 2:s
        [rho_n, rinfo_n] = rt_dm_reconstruct(k0,M0,n0,r);
        if rinfo_n.pval < rinfo.pval
            break;
        end
        rho = rho_n;
        rinfo = rinfo_n;
    end
    return;
end

% Params
eps = 1e-10;
N_max = 100000;
alp = 0.5;

% First approximation by pseudo-inverse
[M,n,k] = rt_vectorize_proto(M0,n0,k0);
[rho, B] = rt_dm_pinv(k, M, n, r);
[U,D] = svd(rho);
c = U*sqrt(D);
c = c(:,1:r);
c = c / sqrt(trace(c'*c));

s = size(rho,1);
I = inv(reshape(B'*n, s, s));

% Iterations
for i = 1:N_max
    cp = c;
    c = alp*get_new_value(c, k, B, I, s) + (1-alp)*cp;
    dc = norm(cp-c);
    if (dc < eps)
        break;
    end
end

rho = c*c';
rho = rho / trace(rho);
[rinfo.pval, rinfo.chi2] = rt_dm_significance(rho, M0, n0, k0);
rinfo.iter = i;
rinfo.r = r;
end

function c = get_new_value(c, k, B, I, s)
rho = c*c';
a = k./(B*rho(:));
J = reshape(B'*a, s, s);
c = I*J*c;
end