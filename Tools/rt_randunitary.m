function u = rt_randunitary(dim)
%RANDUNITARY Generates random unitary matrix of dimension dxd

[q, r] = qr(randn(dim) + 1j * randn(dim));
r = diag(r);
u = q * diag(r ./ abs(r));

end

