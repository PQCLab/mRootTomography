function proc = rt_randprocess(dim, varargin)

op.rank = 1;
op.form = 'chi';
for ja = 1:2:length(varargin)
    op.(lower(varargin{ja})) = varargin{ja + 1};
end

u = rt_randunitary(dim * op.rank);
u = u(:, 1:dim);

proc = permute(reshape(u.', dim, dim, []), [2, 1, 3]);
if strcmpi(op.form, 'kraus')
    return;
end

proc = reshape(proc, dim^2, op.rank);
if strcmpi(op.form, 'root')
    return;
end

proc = proc * proc';

end

