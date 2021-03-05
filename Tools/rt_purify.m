function c = rt_purify(dm, r)
% RT_PURIFY Performs the density matrix purification by taking its square root
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
if nargin < 2
    r = rank(dm);
end

[u, v] = eigs(dm, r);
v = abs(diag(v));
v = rt_to_simplex(v);
c = u * diag(sqrt(v));
end

