classdef TestStateReconstruction < matlab.unittest.TestCase
    properties (Constant)
        Nexp = 500
        Nsample = 1e3
    end
    properties (TestParameter)
        Dim = {2, 2, 3, 3, 3}
        Rank = {1, 2, 1, 2, 3};
    end
    methods(Test, ParameterCombination = 'sequential')
        function testAsymptotic(testCase, Dim, Rank)
            proto = rt_proto_measurement('mub', 'dim', Dim);
            nshots = rt_nshots_divide(1, length(proto), 'equal');
            ex = rt_experiment(Dim, 'state', 'asymptotic').set_data('proto', proto, 'nshots', nshots);

            dm_test = rt_randstate(Dim, 'rank', Rank);
            ex.simulate(dm_test);

            dm_rec = rt_dm_reconstruct(ex, 'Rank', Rank);

            Fidelity = rt_fidelity(dm_test, dm_rec);
            testCase.assertTrue(abs(Fidelity - 1) < 1e-6);
        end
        
        function testPolynomial(testCase, Dim, Rank)
            proto = rt_proto_measurement('mub', 'dim', Dim);
            nshots = rt_nshots_divide(TestStateReconstruction.Nsample, length(proto), 'equal');
            ex = rt_experiment(Dim, 'state').set_data('proto', proto, 'nshots', nshots);
            repeat_run(testCase, ex, Rank, TestStateReconstruction.Nexp, 'StateTestPolynomial');
        end
        
        function testBinomial(testCase, Dim, Rank)
            proto = rt_proto_measurement('mub', 'dim', Dim, 'modifier', 'operator');
            nshots = rt_nshots_divide(TestStateReconstruction.Nsample, length(proto), 'equal');
            ex = rt_experiment(Dim, 'state').set_data('proto', proto, 'nshots', nshots);
            repeat_run(testCase, ex, Rank, TestStateReconstruction.Nexp, 'StateTestBinomial');
        end
        
        function testPoisson(testCase, Dim, Rank)
            proto = rt_proto_measurement('mub', 'dim', Dim, 'modifier', 'operator');
            nshots = rt_nshots_divide(TestStateReconstruction.Nsample, length(proto), 'equal');
            ex = rt_experiment(Dim, 'state', 'poisson_unity').set_data('proto', proto, 'nshots', nshots);
            repeat_run(testCase, ex, Rank, TestStateReconstruction.Nexp, 'StateTestPoisson');
        end
    end
end
