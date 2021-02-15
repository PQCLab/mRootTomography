function proto = rt_proto_process(proto_prep, proto_meas)

if ~iscell(proto_prep)
    dims = size(proto_prep);
    proto_prep = reshape(mat2cell(proto_prep, dims(1), dims(2), ones(1, dims(3))), 1, []);
end
lenp = length(proto_prep);

if ~iscell(proto_meas)
    dims = size(proto_meas);
    proto_meas = reshape(mat2cell(proto_meas, dims(1), dims(2), ones(1, dims(3))), 1, []);
end
lenm = length(proto_meas);

proto_prep = repmat(reshape(proto_prep, 1, []), lenm, 1);
proto_meas = repmat(reshape(proto_meas, 1, []), lenp, 1);

proto = [reshape(proto_prep, [], 1), reshape(proto_meas.', [], 1)].';

end

