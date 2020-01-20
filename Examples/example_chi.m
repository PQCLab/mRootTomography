clc;
clear;
rng(1123);

N = 1; % Number of qubits
r = 2; % Process rank
rr = 2; % Reconstruction rank

s = 2^N;
E = kraus_matrices('random', [s r]);
chi = kraus2chi(E);

proto = rt_chi_protocol('cube', 'pauli', N);
nshots = rt_nshots_devide(1e3,length(proto),'equal');
clicks = rt_simulate(chi, proto, nshots);

[chir,rinfo] = rt_chi_reconstruct(clicks,proto,nshots,'Rank',rr);
rt_fidelity(chir, chi)

