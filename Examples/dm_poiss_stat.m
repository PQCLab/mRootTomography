% Experiment conditions
dim = 2;
r_true = 1;
r_rec = 1;
nshots = 1e3;
proto = cellfun(@(pr) pr(:,:,1), rt_proto_measurement('tetra'), 'UniformOutput', false);
proto = cat(3, proto{:});

% Generate state and data
dm_true = rt_randstate(dim, r_true);
clicks = rt_simulate(dm_true, proto, nshots, 'poiss');

% Reconstruct state and compare to expected one (the true one in our case)
dm_expected = dm_true;
dm_rec = rt_dm_reconstruct(clicks, proto, nshots, 'Rank', r_rec, 'StatType', 'poiss', 'Display', 10);
Fidelity = rt_fidelity(dm_rec, dm_expected);
fprintf('Fidelity: %.6f\n', Fidelity);

% Plot infidelity distribution
d = rt_dm_theory(dm_expected, proto, nshots);
[p, df] = rt_gchi2pdf([], d);
figure; hold on; grid on;
plot(df, p, 'LineWidth', 1.5, 'DisplayName', 'Theory');
plot([1,1] * (1 - Fidelity), ylim, 'LineWidth', 1.5, 'DisplayName', 'Reconstruction');
xlabel('$$1-F$$', 'Interpreter', 'latex');
legend('show');