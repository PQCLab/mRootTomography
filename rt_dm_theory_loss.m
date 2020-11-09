function loss = rt_dm_theory_loss(dm, proto, varargin)
%RT_DM_THEORY_LOSS TODO

nshots = rt_nshots_devide(1, length(proto), 'total');
dF = sum(rt_dm_theory(dm, proto, nshots, varargin{:}));
loss = dF * sum(nshots);

end

