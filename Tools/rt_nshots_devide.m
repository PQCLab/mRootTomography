function nshots = rt_nshots_devide(n, m, method)
%RT_NSHOTS_DEVIDE TODO

if length(n) > 1 || floor(n) ~= n
    error('RT:NshotsDivisionTotal', 'Total shots number should be an integer');
end

if isinf(n)
    nshots = inf(1, m);
    return;
end

if nargin < 3
    method = 'total_int';
end

switch method
    case 'total'
        nshots = ones(1, m) * (n / m);
    case 'total_int'
        nshots = ones(1, m) * floor(n/m);
        idx = 1:(n-sum(nshots));
        nshots(idx) = nshots(idx) + 1;
    case 'equal'
        nshots = ones(1, m) * n;
    otherwise
        error('RT:NshotsDivisionMethod', 'Unknown nshots division method');
end

end

