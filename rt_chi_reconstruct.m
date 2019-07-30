function [chi, rinfo] = rt_chi_reconstruct(k, M, n, r, verbose)
%rt_CHI_RECONSTRUCT Performs the reconstruction of chi-matrix, normalized
%to unity. The reconstruction is by monotone proximal gradient descend based on the paper
%https://papers.nips.cc/paper/5728-accelerated-proximal-gradient-methods-for-nonconvex-programming.pdf
%   INPUT:
%   k - number of counts matrix
%       k(j,i) - number of counts for j-th measurement operator of i-th
%       measurement scheme
%   M - measurement operators for chi-matrix,
%       M{i}(:,:,j) - j-th operator that corresponds to the i-th
%       measurement scheme
%   n - number of measurements vector
%       n(i) - number of measurments fot i-th measurement scheme
%   r - recinstruction rank
%       'auto' - automatic search for adequate rank: starting from rank-1,
%       rank is being increased till the significance level start to
%       decrease
%   verbose - display iterations (default: true)
%
%   OUTPUT:
%   chi - chi-matrix (normalized to unity)
%   rinfo - structure with reconstruction results info
%       rinfo.iter - number of iterations
%       rinfo.logL - logarithmic likelihood
%       rinfo.SL - significance level
%       rinfo.r - resulting rank

if length(n) == 1
    n = ones(length(M),1)*n;
end
if nargin < 5
    verbose = true;
end

% Recursive method for automatic rank estimation
s = size(M{1},1);
if (ischar(r) && strcmpi(r, 'auto')) || r == 0
    [chi, rinfo] = rt_chi_reconstruct(k,M,n,1,verbose);
    for r = 2:s
        [chi_n, rinfo_n] = rt_chi_reconstruct(k,M,n,r,verbose);
        if rinfo_n.SL < rinfo.SL
            break;
        end
        chi = chi_n;
        rinfo = rinfo_n;
    end
    return;
end

% Params
dispfreq = 50 * double(verbose);
eps = 1e-10;
N_max = 100000;
L = sum(n); % Lipschitz constant (tune if there is no convergence)

if dispfreq > 0
    h = rt_fprintreplace('Initialization...');
end

% First approximation by pseudo-inverse
[chi, B, kv] = rt_chi_pinv(k, M, n, r);
[U,S,~] = svd(chi);
e = U*sqrt(S);
e = e(:,1:r);
e = e/sqrt(trace(e'*e));

% We need B, independent of number of measurement
nb = repmat(n,size(k,1),1);
B = B ./ repmat(nb(:),1,s^2);

% Protocol unity
I = 0;
for j = 1:length(M)
    I = I + n(j)*sum(M{j},3);
end

% Accelerated projective gradient descend
x0 = e;
x1 = e;
z1 = e;
t0 = 0;
t1 = 1;
alpx = 1/L;
alpy = 1/L;
F_x1 = 0;
for i = 1:N_max
    y1 = x1 + t0/t1*(z1-x1) + (t0-1)/t1*(x1-x0);
    z2 = project_to_chi(y1-alpy*get_grad(y1, kv, B, I, s), s, r);
    v2 = project_to_chi(x1-alpx*get_grad(x1, kv, B, I, s), s, r);
    F_z2 = get_func(z2,kv,B);
    F_v2 = get_func(v2,kv,B);
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
        h = rt_fprintreplace(sprintf('r = %d | dx = %.4e | dF = %.4e | F = %.4e | %d', r, dx, dF, F_x1, i), h);
    end
    if stopIter
        break;
    end
end
if dispfreq > 0
    fprintf('\n');
end
e = x1/sqrt(trace(x1'*x1));
chi = e*e';

% Reconstruction info
rinfo.iter = i;
rinfo.logL = -get_func(e,kv,B);
[rinfo.SL, rinfo.chi2] = rt_chi_significance(chi,M,n,k);
rinfo.r = r;

end

function F = get_func(e, k, B)
F = -sum(k.*log(get_p(e,B)));
end

function D = get_grad(e, k, B, I, s)
a = k./get_p(e,B);
J = reshape(B'*a, s, s);
D = 2*(I-J)*e;
end

function p = get_p(e, B)
chi = e*e';
chi = chi / trace(chi) * 2;
p = real(B*chi(:));
end

function e = project_to_chi(e, s2, r)

s = sqrt(s2);

Ue = zeros(s*r, s);
for i = 1:r
    Ue((1:s)+(i-1)*s, 1:s) = reshape(e(:,i), s, s)*sqrt(s);
end

[U,S,V] = svd(Ue);
S(S > 0) = 1;
Ue = U*S*V';

for i = 1:r
    E = Ue(((i-1)*s+1):(i*s),1:s);
    e(:,i) = reshape(E, s2, 1)/sqrt(s);
end

end
