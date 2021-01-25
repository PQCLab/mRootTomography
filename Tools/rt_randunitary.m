function u = rt_randunitary(d)
%RANDUNITARY Generates random unitary matrix of dimension sxs

if nargin < 1
    d = 2;
end

[q,r] = qr(randn(d)+1j*randn(d));
r = diag(r);
u = q*diag(r./abs(r));

end

