function [f, x] = rt_gchi2cdf(x, d)
% RT_GCHI2CDF Calculates the generalized chi-squared distribution cumulative distribution function
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
[p, xq] = rt_gchi2pdf([], d);
f = cumsum(p * (xq(2) - xq(1)));
if isempty(x)
    x = xq;
else
    f = interp1(xq, f, x, 'linear', 'extrap');
end
end
