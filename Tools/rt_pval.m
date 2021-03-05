function pval = rt_pval(chi2, df)
% RT_PVAL Calculates the chi-squared test p-value for a specific number of degrees of freedom
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
if df <= 0
    pval = nan;
else
    pval = gammainc(chi2/2,df/2,'upper'); % same as 1-chi2cdf(chi2, df);
    if pval == 0
        pval = vpa(gammainc(chi2/2,df/2,'upper'), 100); % need more precision
    end
end
end

