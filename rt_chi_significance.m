function [SL, chi2] = rt_chi_significance(chi, M, n, k)
%rt_CHI_SIGNIFICANCE Chi-squared test for chi-matrix
%   INPUT:
%   chi - chi-matrix (normalized to unity)
%   M - measurement operators (see description of rt_chi_reconstruct)
%   n - number of measurements (see description of rt_chi_reconstruct)
%   k - number of counts (see description of rt_chi_reconstruct)
%   OUTPUT:
%   SL - significance level
%   chi2 - value of chi-square

mb = numel(k);
mc = length(M);
if mb == mc
    mc = 0;
end

s2 = size(chi, 1);
r = rank(chi);

nE = real(rt_chi_B(M,n)*chi(:));
nO = k(:);

nuP = (2*s2-r)*r - s2;
nu = mb - mc - nuP;

chi2 = sum((nE-nO).^2./nE);
% SL = 1-chi2cdf(chi2, nu);
SL = vpa(gammainc(chi2/2,nu/2,'upper'), 1000);

end
