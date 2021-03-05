function B = rt_meas_matrix(P)
% RT_MEAS_MATRIX Generates the measurement matrix from the measurements operators
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
B = reshape(permute(P,[3,2,1]), size(P,3), []);
end

