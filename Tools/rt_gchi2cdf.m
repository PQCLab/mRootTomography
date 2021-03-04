function [f, x] = rt_gchi2cdf(x, d)
%GCHI2PDF TODO

[p, xq] = rt_gchi2pdf([], d);
f = cumsum(p * (xq(2) - xq(1)));
if isempty(x)
    x = xq;
else
    f = interp1(xq, f, x, 'linear', 'extrap');
end

end

