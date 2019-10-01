function nshots = rt_nshots_devide(n, m, method)
%RT_NSHOTS_DEVIDE TODO

if nargin < 3
    method = 'total_int';
end

switch method
    case 'total'
        nshots = ones(1,m)*(n/m);
    case 'total_int'
        nshots = floor(ones(1,m)*(n/m));
        nshots(end) = n - sum(nshots(1:(end-1)));
    case 'equal'
        nshots = ones(1,m)*n;
    otherwise
        error('Unknown method');
end

end

