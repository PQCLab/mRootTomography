function proc = rt_randprocess(dim, varargin)
% RT_RANDPROCESS Generates a fixed rank quantum process using the extended dynamics representation
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
op.rank = 1;
op.form = 'chi';
op.tracepreserving = true;
for ja = 1:2:length(varargin)
    op.(lower(varargin{ja})) = varargin{ja + 1};
end

if op.tracepreserving
    u = rt_randunitary(dim * op.rank);
    u = u(:, 1:dim);
    proc = permute(reshape(u.', dim, dim, []), [2, 1, 3]);
    switch op.form
        case 'root'
            proc = rt_process_reform(proc, 'kraus2root');
        case 'kraus'
            return;
        case 'chi'
            proc = rt_process_reform(proc, 'kraus2chi');
        otherwise
            error('RT:ProcessForm', 'Unknown process form. Only `chi`, `kraus` and `root` are available.');
    end
else
    proc = randn(dim^2, op.rank) + 1j * randn(dim^2, op.rank);
    proc = proc / sqrt(trace(proc'*proc) / dim);
    switch op.form
        case 'root'
            return;
        case 'kraus'
            proc = rt_process_reform(proc, 'root2kraus');
        case 'chi'
            proc = rt_process_reform(proc, 'root2chi');
        otherwise
            error('RT:ProcessForm', 'Unknown process form. Only `chi`, `kraus` and `root` are available.');
    end
end

end

