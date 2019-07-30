clc;
clear;

N = 1;
s = 2^N;
r = 1;
rho = randstate(s, 'mixed', r);

n = 1e2;
M = rt_protocol_measurement('pauli', N);
[k,n] = rt_dm_simulate(rho, M, n);

[rhor, rinfo] = rt_dm_reconstruct(k,M,n,r);
fidelityB(rho, rhor)