function A = rt_kron3d(A1, A2)
%RT_KRON3D Summary of this function goes here TODO
%   Detailed explanation goes here TODO

A = [];
for i1 = 1:size(A1,3)
    for i2 = 1:size(A2,3)
        A = cat(3, A, kron(A1(:,:,i1),A2(:,:,i2)));
    end
end

end

