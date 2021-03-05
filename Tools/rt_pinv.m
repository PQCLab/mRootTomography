function [dm, c] = rt_pinv(P, freq, r)
% RT_PINV Performs the quantum state reconstruction using Moore-Penrose pseudo-inversion
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
dim = size(P, 1);
if nargin < 3
    r = dim;
end

B = rt_meas_matrix(P);
dm = reshape(mldivide(B, freq), dim, dim);
c = rt_purify(dm, r);
dm = c*c';
end

