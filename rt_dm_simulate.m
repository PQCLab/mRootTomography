function [k,n] = rt_dm_simulate(rho, M, n, asymp)
%rt_SIMULATE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    asymp = false;
end

m = length(M);
s = size(rho, 1);

if length(n) == 1
    n = ones(1,m) * n;
end

k = zeros(s,m);
for j = 1:m
    p = zeros(1,s);
    for i = 1:s
        p(i) = abs(trace(M{j}(:,:,i) * rho));
    end
    if asymp
        k(:,j) = n(j)*p;
    else
        k(:,j) = mnrnd(n(j), p/sum(p))';
    end
end

end

