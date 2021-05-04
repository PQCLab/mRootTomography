% This script performs 500 tomographic numerical experiments with a single
% randomly chosen input state. In each experiment the data is simulated and
% the quantum state is reconstructed using the prior knowledge of the true
% density matrix rank.

% Experiment conditions
n_exp = 500;
dim = 2;
r_true = 1;
r_rec = 1;
nshots = 1e4;
proto = rt_proto_measurement('mub', 'dim', dim);
% proto = rt_proto_measurement('mub', 'dim', dim, 'modifier', 'operator'); % Uncomment to test binomial stats
ex = rt_experiment(dim, 'state').set_data('proto', proto, 'nshots', nshots);

% Generate state
dm_true = rt_randstate(dim, 'Rank', r_true);

% Conduct experiments
Fidelity = zeros(n_exp, 1);
Pval = zeros(n_exp, 1);
for je = 1:n_exp
    fprintf('Experiment %d/%d\n', je, n_exp);
    ex.simulate(dm_true);
    
    [dm_rec, rinfo] = rt_dm_reconstruct(ex, 'Rank', r_rec, 'GetStats', true);
    Fidelity(je) = rt_fidelity(dm_rec, dm_true);
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
d = rt_bound(dm_true, ex);
[p, df] = rt_gchi2pdf([], d);
plot(df, p, 'LineWidth', 1.5, 'DisplayName', 'Theory');
xlabel('$$1-F$$', 'Interpreter', 'latex');
ylabel('$$P(1-F)$$', 'Interpreter', 'latex');
legend('show');