clc;
clear;

N = 2;
s = 2^N;
r = 4;
E = kraus_matrices('random', [s,r]);
chi = kraus2chi(E);

N_exp = 250;
n = ones(1,(6*3)^N)*2048;
[M, P0, M0] = rt_chi_protocol('cube', 'pauli', N);

F = zeros(N_exp, 1);
SL = zeros(N_exp, 1);
for i = 1:N_exp
    disp(i);
    
    k = rt_chi_simulate(chi, P0, M0, n);
    chir = rt_chi_reconstruct(k,M,n,r);
    
    F(i) = fidelityB(chir, chi);
    SL(i) = rt_chi_significance(chir,M,n,k);
end

%%
figure;
histogram(SL);

%%
figure;
h = histogram(1-F);
hold on;
xmax = h.BinLimits(2);
dx = xmax / 100;
x = 0:dx:xmax;
d = rt_chi_theory(chi,M,n);
P = chi2pdf_general(x,d);
N_dist = P * h.BinWidth * N_exp;
plot(x, N_dist, 'LineWidth', 2);

disp([length(d), (2*s^2-r)*r-1, (2*s^2-r)*r-s^2]);
disp([sum(d), mean(1-F)]);
