function l = rt_lossfun(dm, proto, varargin)
% RT_LOSSFUN Calculates the tomography loss function for a quantum state or a quantum process
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
nshots = rt_nshots_divide(1, length(proto), 'total');
dF = sum(rt_bound(dm, proto, nshots, varargin{:}));
l = dF * sum(nshots);
end

