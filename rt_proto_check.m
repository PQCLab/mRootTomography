function [proto,nshots] = rt_proto_check(proto, nshots, clicks)
%RT_NSHOTS_CHECK TODO

if ~iscell(proto)
    proto = num2cell(proto,[1 2]);
end
if nargin > 2 && ischar(nshots) && strcmpi(nshots, 'sum')
    nshots = cellfun(@sum, clicks);
end

if length(nshots) ~= length(proto)
    error('Protocol array should have the same size as nshots array');
end

end

