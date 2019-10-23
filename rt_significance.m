function [pval, chi2, df] = rt_significance(rho, proto, nshots, clicks, r, isChi)
%rt_SIGNIFICANCE TODO

[proto,nshots] = rt_proto_check(proto,nshots);
[M, n, k] = rt_data_join(proto, nshots, clicks);
nE = real(rt_meas_matrix(M)*rho(:)).*n;
nO = k;

s = size(rho,1);
if nargin < 6 || ~isChi
    nuP = (2*s-r)*r - 1;
else
    nuP = (2*s-r)*r - s;
end

sk = cellfun(@(k) length(k), clicks);
df = length(k);
df = df - nuP; % Minus number of the independent real paramenters of the density matrix with rank-r
df = df - sum(sk > 1); % Each measurement with more then 1 output have one connection
df = df - (sum(sk == 1) > 0); % If there exists one measurement with single input there is one more connection (because of the likelihod equation form)

chi2 = sum((nE-nO).^2./nE);

if df == 0
    pval = nan;
else
    pval = gammainc(chi2/2,df/2,'upper'); % same as 1-chi2cdf(chi2, df);
    if pval == 0
        pval = vpa(gammainc(chi2/2,df/2,'upper'), 100); % need more precision
    end
end

end
