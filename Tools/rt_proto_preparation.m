function proto = rt_proto_preparation(ptype, varargin)
%RT_PROTO_PREPARATION Generates preparation protocol
%   INPUT:
%   ptype - protocol type
%       'mub' - Mutually-unbiased bases (MUB) vectors
%       'tetra' - Tetrahedron qubit vectors
%       'octa' - Octahedron qubit vectors
%   varargin - protocol params
%
%   OUTPUT:
%   proto - preparation density matrices
%       proto(:,:,j) - j-th preparation density matrix

switch ptype
    case 'mub'
        dim = varargin{1};
        proto = rt_proto_measurement('mub', dim);
        proto = cat(3, proto{:});
    case 'tetra'
        proto = rt_proto_measurement('tetra', 'operator+');
    case 'octa'
        proto = rt_proto_measurement('tetra', 'operator+-');
    otherwise
        error('RT:UnknownProto', 'Unknown preparation protocol type `%s`', ptype);
end

end

