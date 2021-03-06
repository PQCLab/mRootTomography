function nshots = rt_nshots_divide(n, m, method)
% RT_NSHOTS_DIVIDE Divides a total integer sample size equally over a set of measurements
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
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

