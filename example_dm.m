clc;
clear;

N = 1;
s = 2^N;
r = 1;
rho = rt_randstate(s, 'mixed', r);

nshots = 1e2;
proto = rt_protocol_measurement('pauli', N);
clicks = rt_dm_simulate(rho, proto, nshots);

[rhor, rinfo] = rt_dm_reconstruct(clicks,proto,nshots,r);
rt_fidelity(rho, rhor)