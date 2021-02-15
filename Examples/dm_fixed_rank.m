% Experiment conditions
dim = 2;
r_true = 1;
r_rec = 1;
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
dm_rec = rt_dm_reconstruct(dim, clicks, proto, nshots, 'Rank', r_rec, 'Display', true);
Fidelity = rt_fidelity(dm_rec, dm_expected);
fprintf('Fidelity: %.6f\n', Fidelity);

% Plot infidelity distribution
d = rt_bound(dm_expected, proto, nshots, 'state');
[p, df] = rt_gchi2pdf([], d);
figure; hold on; grid on;
plot(df, p, 'LineWidth', 1.5, 'DisplayName', 'Theory');
plot([1,1] * (1 - Fidelity), ylim, 'LineWidth', 1.5, 'DisplayName', 'Reconstruction');
xlabel('$$1-F$$', 'Interpreter', 'latex');
legend('show');