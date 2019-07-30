clc;
clear;

N = 2;
s = 2^N;
r = 1;
% rho = randstate(s, 'mixed', r);
% rho = genstate('w', 2, true);
% rho = genstate('ghz', 4, true);
c = [1;0;0;exp(1j*pi/4)]/sqrt(2); rho = c*c';

N_exp = 500;
M = rt_protocol_measurement('tetra', N);
n = ones(1,length(M))*1e3;

F = zeros(N_exp, 1);
SL = zeros(N_exp, 1);
Iter = zeros(N_exp, 1);
for i = 1:N_exp
    disp(i);
    
    k = rt_dm_simulate(rho, M, n);
    [rhor, rinfo] = rt_dm_reconstruct(k,M,n,r);
    Iter(i) = rinfo.iter;
    SL(i) = double(rinfo.pval);
    F(i) = fidelityB(rhor, rho);
end

%%
figure;
histogram(Iter);

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
d = rt_dm_theory(rho,M,n);
P = chi2pdf_general(x,d);
N_dist = P * h.BinWidth * N_exp;
plot(x, N_dist, 'LineWidth', 2);
disp([sum(d), mean(1-F)]);
