clc;
clear;
rng(1123);

N = 1; % Number of qubits
r = 2; % Process rank
rr = 2; % Reconstruction rank

s = 2^N;
E = rt_kraus_matrices('random', [s r]);
chi = rt_kraus2chi(E);

n = ones(1,(6*3)^N)*1e3;
[M, P0, M0] = rt_chi_protocol('cube', 'pauli', N);
k = rt_chi_simulate(chi, P0, M0, n);

Lam = [];
for i = 1:length(M)
    for j = 1:size(M{i},3)
        Lam = cat(3,Lam,M{i}(:,:,j));
    end
end
B = rt_chi_B(M,ones(1,(6*3)^N));

[chir,stats] = rt_chi_reconstruct(k,M,n,2);
rt_fidelity(chir, chi)
rt_chi_significance(chir, M, n, k)

