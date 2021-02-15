function c = rt_purify(dm, r)

if nargin < 2
    r = rank(dm);
end

[u, v] = eigs(dm, r);
v = abs(diag(v));
v = v / sum(v);
v = rt_to_simplex(v);
c = u * diag(sqrt(v));

end

