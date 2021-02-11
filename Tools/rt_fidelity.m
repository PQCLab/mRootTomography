function f = rt_fidelity(dm1, dm2)
%RT_FIDELITY Calculates the fidelity of quantum states
%   Calculates the fidelity of quantum states

dm1 = dm1 / trace(dm1);
r1 = rank(dm1);

dm2 = dm2 / trace(dm2);
r2 = rank(dm2);

if r1 == 1 || r2 == 1
    f = abs(trace(dm1 * dm2));
else
    [u, w] = svds(dm1, r1);
    dm1sr = u * sqrt(w) * u';
    a = dm1sr * dm2 * dm1sr;
    f = sum(sqrt(svd(a)))^2;
end

end

