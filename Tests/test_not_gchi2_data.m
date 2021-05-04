function h = test_not_gchi2_data(x, d)

pd = @(x) rt_gchi2cdf(x, d);
h = chi2gof(x(:), 'CDF', pd);
h = logical(h);

end

