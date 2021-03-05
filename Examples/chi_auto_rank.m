% Experiment conditions
dim = 2;
r_true = 1;
nshots = 1e3;
proto_prep = rt_proto_preparation('tetra');
proto_meas = rt_proto_measurement('mub', 'dim', dim);
% proto_meas = rt_proto_measurement('tetra', 'modifier', 'operator'); % Uncomment to test Poisson stats
proto = rt_proto_process(proto_prep, proto_meas);

% Generate state and data
chi_true = rt_randprocess(dim, 'Rank', r_true);
clicks = rt_experiment(dim, 'process')...
    .set_data('proto', proto, 'nshots', nshots)...
    .simulate(chi_true);

% Reconstruct process and compare to the true one
chi_rec = rt_chi_reconstruct(dim, clicks, proto, nshots, 'significanceLevel', 0.01, 'Display', 10);
Fidelity = rt_fidelity(chi_rec, chi_true);
fprintf('Fidelity: %.6f\n', Fidelity);

% Calculate fiducial fidelity bound
d = rt_bound(chi_rec, proto, nshots, 'process');
Fidelity95 = 1 - rt_gchi2inv(0.95, d);
fprintf('Fiducial 95%% fidelity bound: %.6f\n', Fidelity95);

% Plot infidelity distribution
d = rt_bound(chi_true, proto, nshots, 'process');
[p, df] = rt_gchi2pdf([], d);
figure; hold on; grid on;
plot(df, p, 'LineWidth', 1.5, 'DisplayName', 'Theory');
plot([1,1] * (1 - Fidelity), ylim, 'LineWidth', 1.5, 'DisplayName', 'Reconstruction');
xlabel('$$1-F$$', 'Interpreter', 'latex');
legend('show');