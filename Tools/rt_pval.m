function pval = rt_pval(chi2, df)
%RT_PVAL Summary of this function goes here
%   Detailed explanation goes here

if df <= 0
    pval = nan;
else
    pval = gammainc(chi2/2,df/2,'upper'); % same as 1-chi2cdf(chi2, df);
    if pval == 0
        pval = vpa(gammainc(chi2/2,df/2,'upper'), 100); % need more precision
    end
end

end

