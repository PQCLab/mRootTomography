N = 1;
s = 2^N;
r = 2;
% dm_expected = rt_randstate(s, 'mixed', r);
% dm_expected = [1, 1j; -1j, 1]/2;
dm_expected = [0.5, 0.45; 0.45, 0.5];

proto = rt_proto_measurement('pauli', N);
nshots = rt_nshots_devide(1e2,length(proto),'equal');
clicks = rt_simulate(dm_expected, proto, nshots);

[dm, rinfo] = rt_dm_reconstruct(clicks,proto,nshots,'Display',true,'Rank',r);
F = rt_fidelity(dm, dm_expected);
fprintf('Fidelity: %4f\n', F);

d = rt_dm_theory(dm_expected,proto,nshots);
[p,dF] = rt_chi2pdf_general([],d);
figure; hold on; grid on;
plot(dF,p,'LineWidth',1.5,'DisplayName','Theory');
plot([1-F,1-F],ylim,'LineWidth',1.5,'DisplayName','Reconstruction');
xlabel('$$1-F$$','Interpreter','latex');
legend('show');