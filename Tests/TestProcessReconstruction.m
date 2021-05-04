classdef TestProcessReconstruction < matlab.unittest.TestCase
    properties (Constant)
        Nexp = 500
        Nsample = 1e4
    end
    properties (TestParameter)
        Dim = {2, 2, 2, 2, 3, 3, 3, 3, 3}
        Rank = {1, 2, 3, 4, 1, 2, 3, 4, 5};
    end
    methods(Test, ParameterCombination = 'sequential')
        function testAsymptotic(testCase, Dim, Rank)
            proto_prep = rt_proto_preparation('mub', 'dim', Dim);
            proto_meas = rt_proto_measurement('mub', 'dim', Dim);
            proto = rt_proto_process(proto_prep, proto_meas);
            nshots = rt_nshots_divide(1, size(proto, 2), 'equal');
            ex = rt_experiment(Dim, 'process', 'asymptotic').set_data('proto', proto, 'nshots', nshots);

            chi_test = rt_randprocess(Dim, 'rank', Rank);
            ex.simulate(chi_test);

            chi_rec = rt_chi_reconstruct(ex, 'Rank', Rank);

            Fidelity = rt_fidelity(chi_test, chi_rec);
            testCase.assertTrue(abs(Fidelity - 1) < 1e-6);
        end
        
        function testPolynomial(testCase, Dim, Rank)
            proto_prep = rt_proto_preparation('mub', 'dim', Dim);
            proto_meas = rt_proto_measurement('mub', 'dim', Dim);
            proto = rt_proto_process(proto_prep, proto_meas);
            nshots = rt_nshots_divide(TestStateReconstruction.Nsample, size(proto, 2), 'equal');
            ex = rt_experiment(Dim, 'process').set_data('proto', proto, 'nshots', nshots);
            repeat_run(testCase, ex, Rank, TestStateReconstruction.Nexp, 'ProcessTestPolynomial');
        end
        
        function testBinomial(testCase, Dim, Rank)
            proto_prep = rt_proto_preparation('mub', 'dim', Dim);
            proto_meas = rt_proto_measurement('mub', 'dim', Dim, 'modifier', 'operator');
            proto = rt_proto_process(proto_prep, proto_meas);
            nshots = rt_nshots_divide(TestStateReconstruction.Nsample, size(proto, 2), 'equal');
            ex = rt_experiment(Dim, 'process').set_data('proto', proto, 'nshots', nshots);
            repeat_run(testCase, ex, Rank, TestStateReconstruction.Nexp, 'ProcessTestBinomial');
        end
        
        function testPoisson(testCase, Dim, Rank)
            proto_prep = rt_proto_preparation('mub', 'dim', Dim);
            proto_meas = rt_proto_measurement('mub', 'dim', Dim, 'modifier', 'operator');
            proto = rt_proto_process(proto_prep, proto_meas);
            nshots = rt_nshots_divide(TestStateReconstruction.Nsample, size(proto, 2), 'equal');
            ex = rt_experiment(Dim, 'process', 'poisson').set_data('proto', proto, 'nshots', nshots);
            repeat_run(testCase, ex, Rank, TestStateReconstruction.Nexp, 'ProcessTestPoisson');
        end
    end
end
