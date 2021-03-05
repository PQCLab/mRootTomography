function dm = rt_randstate(dim, varargin)
% RT_RANDSTATE Generates a fixed rank quantum state using the partial tracing
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
op.rank = dim;
for ja = 1:2:length(varargin)
    op.(lower(varargin{ja})) = varargin{ja + 1};
end

c = randn(dim, op.rank) + 1j * randn(dim, op.rank);
dm = c * c';
dm = dm / trace(dm);
end

