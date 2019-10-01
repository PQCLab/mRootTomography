function rt_proto_check(proto, nshots)
%RT_NSHOTS_CHECK TODO

if length(nshots) ~= length(proto)
    error('Protocol array should have the same size as nshots array');
end

end

