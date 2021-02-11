% Experiment conditions
dim = 4;
r_true = 2;
nshots = 1e3;
proto = rt_proto_measurement('mub', dim);

% Generate state and data
dm_true = rt_randstate(dim, r_true);
clicks = rt_simulate(dm_true, proto, nshots);

% Reconstruct state and compare to expected one (the true one in our case)
dm_expected = dm_true;
[dm_rec, rinfo] = rt_dm_reconstruct(clicks, proto, nshots, 'significanceLevel', 0.01, 'Display', 10);
Fidelity = rt_fidelity(dm_rec, dm_expected);
fprintf('Fidelity: %.6f\n', Fidelity);