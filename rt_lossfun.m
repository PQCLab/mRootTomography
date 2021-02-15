function loss = rt_lossfun(dm, proto, varargin)
%RT_LOSS TODO

nshots = rt_nshots_devide(1, length(proto), 'total');
dF = sum(rt_bound(dm, proto, nshots, varargin{:}));
loss = dF * sum(nshots);

end

