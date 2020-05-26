function [proto,nshots] = rt_proto_check(proto, nshots, clicks)
%RT_NSHOTS_CHECK TODO

if ~iscell(proto)
    proto = num2cell(proto,[1 2]);
end
if nargin > 2 && ischar(nshots) && strcmpi(nshots, 'sum')
    if ischar(nshots) && strcmpi(nshots, 'sum')
        lens = cellfun(@length, clicks);
        if sum(lens == 1) > 0
            error('You must specify number of experiment repetitions for experiments with a single possible outcome');
        end
        nshots = cellfun(@sum, clicks);
    end
    if length(nshots) ~= length(proto)
        error('Protocol array should have the same size as nshots array');
    end
end

end

