function ne = rt_chi_expval(chi, M, n)
%rt_CHI_EXPVAL Summary of this function goes here
%   Detailed explanation goes here

m = length(M);
s = size(M{1},1);
if m ~= length(n)
    error('Wrong number of measurements');
end

ne = zeros(s,m);
for j = 1:m
    for i = 1:s
        ne(i,j) = real(trace(M{j}(:,:,i)*chi))*s*n(j);
    end
end

end

