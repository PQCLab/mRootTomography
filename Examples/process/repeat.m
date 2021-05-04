% This script performs 500 tomographic numerical experiments with a single
% randomly chosen process. In each experiment the data is simulated and
% the quantum process matrix is reconstructed using the prior knowledge
% of the true process matrix.

% Experiment conditions
n_exp = 500;
dim = 2;
r_true = 1;
r_rec = 1;
nshots = 1e4;
proto_prep = rt_proto_preparation('mub', 'dim', dim);
proto_meas = rt_proto_measurement('mub', 'dim', dim);
% proto_meas = rt_proto_measurement('mub', 'dim', dim, 'modifier', 'operator'); % Uncomment to test binomial stats
proto = rt_proto_process(proto_prep, proto_meas);
ex = rt_experiment(dim, 'process').set_data('proto', proto, 'nshots', nshots);

% Generate state
chi_true = rt_randprocess(dim, 'Rank', r_true);

% Conduct experiments
Fidelity = zeros(n_exp,1);
Pval = zeros(n_exp,1);
for je = 1:n_exp
    fprintf('Experiment %d/%d\n', je, n_exp);
    ex.simulate(chi_true);
    
    [chi_rec, rinfo] = rt_chi_reconstruct(ex, 'Rank', r_rec, 'getStats', true);
    Fidelity(je) = rt_fidelity(chi_rec, chi_true);
    Pval(je) = double(rinfo.pval);
end

%% Plot results
figure;
hold on; grid on;
histogram(Pval, 'Normalization', 'pdf', 'DisplayName', 'Numerical Experiments');
plot([0, 1], [1, 1], 'LineWidth', 1.5, 'DisplayName', 'Theory');
xlabel('$$p-value$$', 'Interpreter', 'latex');
ylabel('$$P(p-value)$$', 'Interpreter', 'latex');
legend('show');

figure;
hold on; grid on;
histogram(1 - Fidelity, 'Normalization', 'pdf', 'DisplayName', 'Numerical Experiments');
d = rt_bound(chi_true, ex);
[p, df] = rt_gchi2pdf([], d);
plot(df, p, 'LineWidth', 1.5, 'DisplayName', 'Theory');
xlabel('$$1-F$$', 'Interpreter', 'latex');
ylabel('$$P(1-F)$$', 'Interpreter', 'latex');
legend('show');