function dm = rt_randstate(dim, varargin)
%RANDSTATE Generates random state
%   Generates random d-dimention density matrix of rank r
op.rank = dim;
for ja = 1:2:length(varargin)
    op.(lower(varargin{ja})) = varargin{ja + 1};
end

c = randn(dim, op.rank) + 1j * randn(dim, op.rank);
dm = c * c';
dm = dm / trace(dm);

end

