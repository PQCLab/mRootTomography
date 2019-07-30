function [M,n,k] = rt_vectorize_proto(M,n,k)

if ~iscell(M)
    return;
end

M0 = M;
M = [];
for j = 1:length(M0)
    M = cat(3,M,M0{j});
end
n = reshape(repmat(n, size(k,1), 1),[],1);
k = k(:);

end

