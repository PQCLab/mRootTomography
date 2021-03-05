function A = rt_kron3d(A1, A2)
% RT_KRON_3D Calculates the generalized Kronecker product of two 3D arrays
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
A = [];
for i1 = 1:size(A1,3)
    for i2 = 1:size(A2,3)
        A = cat(3, A, kron(A1(:,:,i1),A2(:,:,i2)));
    end
end
end

