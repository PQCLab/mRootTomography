function [dm, psi] = rt_pinv(x, freq, r)
% RT_PINV Performs the quantum state reconstruction using Moore-Penrose pseudo-inversion
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
if size(x, 3) == 1
    bmat = x;
else
    bmat = rt_meas_matrix(x);
end

dim = fix(sqrt(size(bmat, 2)));
if nargin < 3
    r = dim;
end

dm = reshape(mldivide(bmat, freq), dim, dim);
psi = rt_purify(dm, r);
dm = psi * psi';
end

