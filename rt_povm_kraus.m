function M = rt_povm_kraus(M0, E)

M = zeros(size(M0));
for j = 1:size(M0,3)
    for k = 1:size(E,3)
        M(:,:,j) = M(:,:,j) + E(:,:,k)'*M0(:,:,j)*E(:,:,k);
    end
end

end

