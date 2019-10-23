clc;
clear;

N = 3;
s = 2^N;
r = 5;
rho = rt_randstate(s, 'mixed', r);

N_exp = 100;
proto = rt_proto_measurement('pauli', N);
nshots = rt_nshots_devide(1e5, length(proto));

F = zeros(N_exp,1);
Pval = zeros(N_exp,1);
for i = 1:N_exp
    fprintf('Experiment %d/%d\n', i, N_exp);
    
    clicks = rt_simulate(rho, proto, nshots);
    [rhor, rinfo] = rt_dm_reconstruct(clicks,proto,nshots,'Rank',r);
    Pval(i) = double(rinfo.pval);
    F(i) = rt_fidelity(rhor, rho);
end

%%
figure;
histogram(Pval);

%%
figure;
h = histogram(1-F);
hold on;
xmax = h.BinLimits(2);
dx = xmax / 100;
x = 0:dx:xmax;
d = rt_dm_theory(rho,proto,nshots,r);
P = rt_chi2pdf_general(x,d);
N_dist = P * h.BinWidth * N_exp;
plot(x, N_dist, 'LineWidth', 2);
disp([sum(d), mean(1-F)]);
