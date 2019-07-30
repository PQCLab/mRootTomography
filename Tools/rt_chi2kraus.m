function E = rt_chi2kraus(chi, tol)
%CHI2KRAUS Calculates kraus matrices from chi-matrix
%   Chi-matrix should be:
%       - normalized to unity
%       - in natural representation

if nargin < 2
    tol = 1e-10;
end

s = sqrt(size(chi,1));
[U,D] = svd(chi/trace(chi)*s);
D = diag(D);
ind = find(D > tol)';

E = zeros(s,s,length(ind));
for i = 1:length(ind)
    E(:,:,i) = reshape(U(:,ind(i))*sqrt(D(ind(i))), s, s);
end

end

