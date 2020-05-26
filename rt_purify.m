function c = rt_purify(dm, r)

[U,D] = eig(dm);
[d,ind] = sort(rt_projprob(real(diag(D))),'descend');
U = U(:,ind);
if nargin > 1
    d = d(1:r);
    U = U(:,1:r);
end
c = U*diag(sqrt(d));

end

