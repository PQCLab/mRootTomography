function l = rt_lossfun(dm, proto, varargin)
%RT_LOSS TODO

nshots = rt_nshots_divide(1, length(proto), 'total');
dF = sum(rt_bound(dm, proto, nshots, varargin{:}));
l = dF * sum(nshots);

end

