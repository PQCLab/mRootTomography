% This script simulates the data and performs the quantum process
% reconstruction using the prioir knowledge of the true process matrix rank

% Experiment conditions
dim = 2;
r_true = 1;
r_rec = 1;
nshots = 1e4;
proto_prep = rt_proto_preparation('mub', 'dim', dim);
% proto_meas = rt_proto_measurement('mub', 'dim', dim);
proto_meas = rt_proto_measurement('mub', 'dim', dim, 'modifier', 'operator'); % Uncomment to test binomial stats
proto = rt_proto_process(proto_prep, proto_meas);
ex = rt_experiment(dim, 'process').set_data('proto', proto, 'nshots', nshots);

% Generate state and data
chi_true = rt_randprocess(dim, 'Rank', r_true);
ex.simulate(chi_true);

% Reconstruct process and compare to expected one (the true one in our case)
chi_rec = rt_chi_reconstruct(ex, 'Rank', r_rec, 'Display', 10);
Fidelity = rt_fidelity(chi_rec, chi_true);
fprintf('Fidelity: %.6f\n', Fidelity);

% Calculate fiducial fidelity bound
d = rt_bound(chi_rec, ex);
Fidelity95 = 1 - rt_gchi2inv(0.95, d);
fprintf('Fiducial 95%% fidelity bound: %.6f\n', Fidelity95);

% Plot infidelity distribution
d = rt_bound(chi_true, ex);
[p, df] = rt_gchi2pdf([], d);
figure; hold on; grid on;
plot(df, p, 'LineWidth', 1.5, 'DisplayName', 'Theory');
plot([1,1] * (1 - Fidelity), ylim, 'LineWidth', 1.5, 'DisplayName', 'Reconstruction');
xlabel('$$1-F$$', 'Interpreter', 'latex');
legend('show');
