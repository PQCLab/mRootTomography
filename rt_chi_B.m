function B = rt_chi_B(M,n)
%rt_PINV Returns measurement matrix (see
%https://arxiv.org/pdf/1106.2056.pdf) for the reconstruction of chi-matrix
%   INPUT:
%   M - measurement operators (see description of rt_chi_reconstruct)
%   n - number of measurements (see description of rt_chi_reconstruct)
%
%   OUTPUT:
%   B - measurement matrix

m = length(n);
mproj = size(M{1},3);
s2 = size(M{1},1);
s = sqrt(s2);

B = zeros(mproj*m, s2^2);
for j = 1:m
    ind = (1:mproj)+(j-1)*mproj;
    B(ind,:) = n(j)*rt_dm_B(M{j})*s;
end

end

