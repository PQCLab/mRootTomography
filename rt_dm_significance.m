function [SL, chi2] = rt_dm_significance(rho, M, n, k)
%rt_CHI_SIGNIFICANCE Chi-squared test for chi-matrix
%   INPUT:
%   rho - density matrix
%   M - measurement operators (see description of rt_dm_reconstruct)
%   n - number of measurements (see description of rt_dm_reconstruct)
%   k - number of counts (see description of rt_dm_reconstruct)
%   OUTPUT:
%   SL - significance level
%   chi2 - value of chi-square

mb = numel(k);
mc = length(M);
if mb == mc
    mc = 0;
end

[M,n,k] = rt_vectorize_proto(M,n,k);

s = size(rho, 1);
r = rank(rho);

nE = real(rt_dm_B(M,n)*rho(:));
nO = k(:);

nuP = (2*s-r)*r - 1;
nu = mb - mc - nuP;

chi2 = sum((nE-nO).^2./nE);
% SL = 1-chi2cdf(chi2, nu);
SL = vpa(gammainc(chi2/2,nu/2,'upper'), 1000);

end
