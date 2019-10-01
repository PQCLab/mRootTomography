function Loss = rt_dm_theory_loss(dm, proto, r)
%RT_DM_THEORY_LOSS Summary of this function goes here TODO
%   Detailed explanation goes here TODO

nshots = rt_dm_protocol_check(proto, 1e6);
dF = sum(rt_dm_theory(dm, proto, nshots, r));
Loss = dF * sum(nshots);

end

