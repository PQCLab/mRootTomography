function x = rt_gchi2inv(f, d)
% RT_GCHI2INV Calculates the generalized chi-squared distribution inverse cumulative distribution function
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
[p, xq] = rt_gchi2pdf([], d);
fq = cumsum(p * (xq(2) - xq(1)));
idx0 = find(fq > 1e-10, 1);
idx1 = find(fq < 1-1e-10, 1, 'last');
xq = xq(idx0:idx1);
fq = fq(idx0:idx1);
[fq, idx] = unique(fq);
xq = xq(idx);
x = interp1(fq, xq, f, 'linear', 'extrap');
end

