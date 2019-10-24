N = 1;
s = 2^N;
r = 2;
dm_true = rt_randstate(s, 'mixed', r);

proto = rt_proto_measurement('pauli', N);
nshots = rt_nshots_devide(1e2,length(proto),'equal');
clicks = rt_simulate(dm_true, proto, nshots);

dm_expected = dm_true;
[dm, rinfo] = rt_dm_reconstruct(clicks,proto,nshots,'Display',true);
F = rt_fidelity(dm, dm_expected);
fprintf('Fidelity: %4f\n', F);

d = rt_dm_theory(dm_expected,proto,nshots);
[p,dF] = rt_chi2pdf_general([],d);
figure; hold on; grid on;
plot(dF,p,'LineWidth',1.5,'DisplayName','Theory');
plot([1-F,1-F],ylim,'LineWidth',1.5,'DisplayName','Reconstruction');
xlabel('$$1-F$$','Interpreter','latex');
legend('show');