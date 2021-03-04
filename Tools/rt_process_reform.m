function proc = rt_process_reform(proc, desc, varargin)

switch desc
    case 'chi2kraus'
        e = rt_process_reform(proc, 'chi2root', varargin{:});
        proc = rt_process_reform(e, 'root2kraus');
    case 'chi2root'
        if isempty(varargin)
            r = rank(proc);
        else
            r = varargin{1};
        end
        proc = rt_purify(proc, r) * sqrt(trace(proc));
    case 'root2chi'
        proc = proc * proc';
    case 'root2kraus'
        [dim2, r] = size(proc);
        dim = sqrt(dim2);
        proc = reshape(proc, [dim, dim, r]);
    case 'kraus2chi'
        e = rt_process_reform(proc, 'kraus2root');
        proc = rt_process_reform(e, 'root2chi');
    case 'kraus2root'
        [dim, ~, r] = size(proc);
        proc = reshape(proc, [dim^2, r]);
    otherwise
        error('RT:ProcessReformDesc', 'Unknown process reform description');
end

end

