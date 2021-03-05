function u = rt_randunitary(dim)
% RT_RANDUNITARY Generates a Haar random unitary matrix
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
[q, r] = qr(randn(dim) + 1j * randn(dim));
r = diag(r);
u = q * diag(r ./ abs(r));
end
