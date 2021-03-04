function x = rt_gchi2inv(f, d)
%GCHI2PDF TODO

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

