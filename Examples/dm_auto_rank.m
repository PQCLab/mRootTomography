% Experiment conditions
dim = 2;
r_true = 1;
nshots = 1e3;
proto = rt_proto_measurement('mub', dim);
% proto = rt_proto_measurement('tetra', 'operator+-'); % Uncomment to test Poisson stats

% Generate state and data
dm_true = rt_randstate(dim, 'Rank', r_true);
clicks = rt_experiment(dim, 'state')...
    .set_data('proto', proto, 'nshots', nshots)...
    .simulate(dm_true);

% Reconstruct state and compare to expected one (the true one in our case)
dm_expected = dm_true;
[dm_rec, rinfo] = rt_dm_reconstruct(dim, clicks, proto, nshots, 'significanceLevel', 0.01, 'Display', 10);
Fidelity = rt_fidelity(dm_rec, dm_expected);
fprintf('Fidelity: %.6f\n', Fidelity);