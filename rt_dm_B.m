function B = rt_dm_B(M,n)

s = size(M,1);
m = size(M,3);

if nargin < 2
    n = ones(m,1);
end

B = zeros(m, s^2);
for j = 1:m
    B(j,:) = n(j)*reshape(M(:,:,j), s^2, 1)';
end

end

