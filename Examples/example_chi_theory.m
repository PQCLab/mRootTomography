clc;
clear;

N = 1;
s = 2^N;
s2 = s^2;
r = 3;

U = rt_randunitary(s*r);
e = reshape(permute(reshape(transpose(U(:,1:s)),[s,s,r]),[2,1,3]),[s2,r]);
chi = e*e';

N_exp = 250;
proto = rt_chi_protocol('cube', 'pauli', N);
nshots = rt_nshots_devide(1e6, length(proto));

F = zeros(N_exp, 1);
Pval = zeros(N_exp, 1);
for i = 1:N_exp
    fprintf('Experiment %d/%d\n', i, N_exp);
    
    clicks = rt_simulate(chi, proto, nshots);
    [chir, rinfo] = rt_chi_reconstruct(clicks,proto,nshots,'Rank',r);
    Pval(i) = rinfo.pval;
    F(i) = rt_fidelity(chir, chi);
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
d = rt_chi_theory(chi,proto,nshots);
P = chi2pdf_general(x,d);
N_dist = P * h.BinWidth * N_exp;
plot(x, N_dist, 'LineWidth', 2);

disp([length(d), (2*s^2-r)*r-1, (2*s^2-r)*r-s^2]);
disp([sum(d), mean(1-F)]);
