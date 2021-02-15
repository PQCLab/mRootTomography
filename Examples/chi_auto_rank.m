% Experiment conditions
n_exp = 500;
dim = 2;
r_true = 3;
nshots = 1e5;
proto_prep = rt_proto_preparation('tetra');
proto_meas = rt_proto_measurement('mub', dim);
% proto_meas = rt_proto_measurement('tetra', 'operator+-'); % Uncomment to test Poisson stats
proto = rt_proto_process(proto_prep, proto_meas);

% Generate state and data
chi_true = rt_randprocess(dim, 'Rank', r_true);
clicks = rt_experiment(dim, 'process')...
    .set_data('proto', proto, 'nshots', nshots)...
    .simulate(chi_true);

% Reconstruct state and compare to expected one (the true one in our case)
chi_expected = chi_true;
[chi_rec, rinfo] = rt_chi_reconstruct(dim, clicks, proto, nshots, 'significanceLevel', 0.01, 'Display', 10);
Fidelity = rt_fidelity(chi_rec, chi_expected);
fprintf('Fidelity: %.6f\n', Fidelity);