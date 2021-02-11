function dm = rt_randstate(d, r)
%RANDSTATE Generates random state
%   Generates random d-dimention density matrix of rank r
if nargin < 2
    r = d;
end

c = randn(d, r) + 1j * randn(d, r);
dm = c * c';
dm = dm / trace(dm);

end

