% This script simulates the data and performs the quantum state
% reconstruction using the prioir knowledge of the true density matrix rank

% Experiment conditions
dim = 2;
r_true = 1;
r_rec = 1;
nshots = 1e4;
proto = rt_proto_measurement('mub', 'dim', dim);
% proto = rt_proto_measurement('mub', 'dim', dim, 'modifier', 'operator'); % Uncomment to test binomial stats
ex = rt_experiment(dim, 'state').set_data('proto', proto, 'nshots', nshots);

% Generate state and data
dm_true = rt_randstate(dim, 'Rank', r_true);
ex.simulate(dm_true);

% Reconstruct state and compare to the true one
dm_rec = rt_dm_reconstruct(ex, 'Rank', r_rec, 'Display', 10);
Fidelity = rt_fidelity(dm_rec, dm_true);
fprintf('Fidelity: %.6f\n', Fidelity);

% Calculate fiducial fidelity bound
d = rt_bound(dm_rec, ex);
Fidelity95 = 1 - rt_gchi2inv(0.95, d);
fprintf('Fiducial 95%% fidelity bound: %.6f\n', Fidelity95);

% Plot infidelity distribution
d = rt_bound(dm_true, ex);
[p, df] = rt_gchi2pdf([], d);
figure; hold on; grid on;
plot(df, p, 'LineWidth', 1.5, 'DisplayName', 'Theory');
plot([1,1] * (1 - Fidelity), ylim, 'LineWidth', 1.5, 'DisplayName', 'Reconstruction');
xlabel('$$1-F$$', 'Interpreter', 'latex');
legend('show');