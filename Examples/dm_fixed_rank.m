% Experiment conditions
dim = 2;
r_true = 1;
r_rec = 1;
nshots = 1e3;
proto = rt_proto_measurement('mub', 'dim', dim);
% proto = rt_proto_measurement('tetra', 'modifier', 'operator'); % Uncomment to test Poisson stats

% Generate state and data
dm_true = rt_randstate(dim, 'Rank', r_true);
clicks = rt_experiment(dim, 'state')...
    .set_data('proto', proto, 'nshots', nshots)...
    .simulate(dm_true);

% Reconstruct state and compare to the true one
dm_rec = rt_dm_reconstruct(dim, clicks, proto, nshots, 'Rank', r_rec, 'Display', true);
Fidelity = rt_fidelity(dm_rec, dm_true);
fprintf('Fidelity: %.6f\n', Fidelity);

% Calculate fiducial fidelity bound
d = rt_bound(dm_rec, proto, nshots, 'state');
Fidelity95 = 1 - rt_gchi2inv(0.95, d);
fprintf('Fiducial 95%% fidelity bound: %.6f\n', Fidelity95);

% Plot infidelity distribution
d = rt_bound(dm_true, proto, nshots, 'state');
[p, df] = rt_gchi2pdf([], d);
figure; hold on; grid on;
plot(df, p, 'LineWidth', 1.5, 'DisplayName', 'Theory');
plot([1,1] * (1 - Fidelity), ylim, 'LineWidth', 1.5, 'DisplayName', 'Reconstruction');
xlabel('$$1-F$$', 'Interpreter', 'latex');
legend('show');