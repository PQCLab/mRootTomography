function h = test_not_uniform_data(x)

pd = makedist('Uniform');
h = chi2gof(x(:), 'CDF', pd);
h = logical(h);

end