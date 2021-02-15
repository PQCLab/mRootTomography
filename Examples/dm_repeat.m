% Experiment conditions
n_exp = 500;
dim = 2;
r_true = 1;
r_rec = 1;
nshots = 1e5;
proto = rt_proto_measurement('mub', dim);
% proto = rt_proto_measurement('tetra', 'operator+-'); % Uncomment to test Poisson stats

% Generate state
dm_true = rt_randstate(dim, 'Rank', r_true);

% Conduct experiments
dm_expected = dm_true;
Fidelity = zeros(n_exp,1);
Pval = zeros(n_exp,1);
for je = 1:n_exp
    fprintf('Experiment %d/%d\n', je, n_exp);
    clicks = rt_experiment(dim, 'state')...
        .set_data('proto', proto, 'nshots', nshots)...
        .simulate(dm_true);
    
    [dm_rec, rinfo] = rt_dm_reconstruct(dim, clicks, proto, nshots, 'Rank', r_rec, 'getStats', true);
    Fidelity(je) = rt_fidelity(dm_rec, dm_expected);
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
d = rt_bound(dm_expected, proto, nshots, 'state');
[p, df] = rt_gchi2pdf([], d);
plot(df, p, 'LineWidth', 1.5, 'DisplayName', 'Theory');
xlabel('$$1-F$$', 'Interpreter', 'latex');
ylabel('$$P(1-F)$$', 'Interpreter', 'latex');
legend('show');