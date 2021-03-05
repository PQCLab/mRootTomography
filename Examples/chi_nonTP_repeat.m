% Experiment conditions
n_exp = 500;
dim = 2;
r_true = 1;
r_rec = 1;
nshots = 1e3;
proto_prep = rt_proto_preparation('tetra');
proto_meas = rt_proto_measurement('mub', 'dim', dim);
proto = rt_proto_process(proto_prep, proto_meas);

% Generate state
chi_true = rt_randprocess(dim, 'Rank', r_true, 'TracePreserving', false);

% Conduct experiments
chi_expected = chi_true;
Fidelity = zeros(n_exp,1);
Pval = zeros(n_exp,1);
for je = 1:n_exp
    fprintf('Experiment %d/%d\n', je, n_exp);
    clicks = rt_experiment(dim, 'process', 'poiss')...
        .set_data('proto', proto, 'nshots', nshots)...
        .simulate(chi_true);
    
    [chi_rec, rinfo] = rt_chi_reconstruct(dim, clicks, proto, nshots, 'Rank', r_rec, 'getStats', true, 'TracePreserving', false, 'StatType', 'poiss');
    Fidelity(je) = rt_fidelity(chi_rec, chi_expected);
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
d = rt_bound(chi_expected, proto, nshots, 'process', 'TracePreserving', false);
[p, df] = rt_gchi2pdf([], d);
plot(df, p, 'LineWidth', 1.5, 'DisplayName', 'Theory');
xlabel('$$1-F$$', 'Interpreter', 'latex');
ylabel('$$P(1-F)$$', 'Interpreter', 'latex');
legend('show');