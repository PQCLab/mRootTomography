function repeat_run(testCase, ex, Rank, n_exp, name)

Dim = ex.dim;
if strcmp(ex.obj_type, 'state')
    x_test = rt_randstate(Dim, 'rank', Rank);
else
    x_test = rt_randprocess(Dim, 'rank', Rank);
end

Fidelity = zeros(1, n_exp);
Pval = zeros(1, n_exp);
for je = 1:n_exp
    fprintf('%s ==> Dim = %d, Rank = %d, Experiment %d/%d\n', name, ex.dim, Rank, je, n_exp);
    ex.simulate(x_test);

    if strcmp(ex.obj_type, 'state')
        [x_rec, rinfo] = rt_dm_reconstruct(ex, 'Rank', Rank, 'getStats', true);
    else
        [x_rec, rinfo] = rt_chi_reconstruct(ex, 'Rank', Rank, 'getStats', true);
    end
    
    
    Fidelity(je) = rt_fidelity(x_rec, x_test);
    Pval(je) = double(rinfo.pval);
end

if strcmpi(name, 'StateTestPolynomial') && Rank == Dim
    testCase.assertTrue(all(isnan(Pval)));
else
    testCase.assertFalse(test_not_uniform_data(Pval));
end
d = rt_bound(x_test, ex);
testCase.assertFalse(test_not_gchi2_data(1 - Fidelity, d));

end