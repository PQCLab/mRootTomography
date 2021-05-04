function proto = rt_proto_preparation(ptype, varargin)
% RT_PROTO_PREPARATION Generates a set of density matrices according to a specific protocol
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
op.dim = nan;
op.num = nan;
op.nsub = 1;
for ja = 1:2:length(varargin)
    op.(lower(varargin{ja})) = varargin{ja + 1};
end

switch ptype
    case 'mub'
        proto = rt_proto_measurement('mub', 'Dim', op.dim, 'Modifier', 'operator');
    case 'tetra'
        proto = rt_proto_measurement('tetra', 'Modifier', 'operator1');
    case 'octa'
        proto = rt_proto_measurement('tetra', 'Modifier', 'operator');
    case 'random'
        proto = rt_proto_measurement('random_projectors', 'Dim', op.dim, 'Num', op.num);
    otherwise
        error('RT:UnknownProto', 'Unknown preparation protocol type `%s`', ptype);
end

proto = cat(3, proto{:});
if op.nsub > 1
    proto0 = proto;
    for jn = 2:op.nsub
        proto = rt_kron3d(proto, proto0);
    end
end

end

