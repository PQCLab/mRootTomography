function proto = rt_proto_measurement(ptype, varargin)
%RT_PROTO_MEASUREMENT Library of quantum states measurements protocols
%   INPUT:
%   ptype - protocol type
%       'mub' - Mutually-unbiased bases (MUB) protocol
%       'tetra' - Tetrahedron qubit bases
%   varargin - protocol params
%
%   OUTPUT:
%   proto - measurement operators
%       proto{j}(:,:,k) - k-th measurement operator corresponding to the
%       j-th measurement

switch ptype
    case 'mub'
        d = varargin{1};
        pname = ['mub', num2str(d)];
        mubs = load('mubs.mat', pname);
        mubs = mubs.(pname);
        proto = cell(1, d + 1);
        for j = 1:(d+1)
            proto{j} = zeros(d, d, d);
            for k = 1:d
                proto{j}(:,:,k) = mubs(:,k,j) * mubs(:,k,j)';
            end
        end
    case 'tetra'
        proto = cell(1,4);
        proto{1}(:,:,1) = [sqrt(3)+1 sqrt(2)*exp(1j*1*pi/4); sqrt(2)*exp(-1j*1*pi/4) sqrt(3)-1]/(2*sqrt(3));
        proto{1}(:,:,2) = [sqrt(3)-1 sqrt(2)*exp(1j*5*pi/4); sqrt(2)*exp(-1j*5*pi/4) sqrt(3)+1]/(2*sqrt(3));
        proto{2}(:,:,1) = [sqrt(3)-1 sqrt(2)*exp(1j*7*pi/4); sqrt(2)*exp(-1j*7*pi/4) sqrt(3)+1]/(2*sqrt(3));
        proto{2}(:,:,2) = [sqrt(3)+1 sqrt(2)*exp(1j*3*pi/4); sqrt(2)*exp(-1j*3*pi/4) sqrt(3)-1]/(2*sqrt(3));
        proto{3}(:,:,1) = [sqrt(3)+1 sqrt(2)*exp(1j*5*pi/4); sqrt(2)*exp(-1j*5*pi/4) sqrt(3)-1]/(2*sqrt(3));
        proto{3}(:,:,2) = [sqrt(3)-1 sqrt(2)*exp(1j*1*pi/4); sqrt(2)*exp(-1j*1*pi/4) sqrt(3)+1]/(2*sqrt(3));
        proto{4}(:,:,1) = [sqrt(3)-1 sqrt(2)*exp(1j*3*pi/4); sqrt(2)*exp(-1j*3*pi/4) sqrt(3)+1]/(2*sqrt(3));
        proto{4}(:,:,2) = [sqrt(3)+1 sqrt(2)*exp(1j*7*pi/4); sqrt(2)*exp(-1j*7*pi/4) sqrt(3)-1]/(2*sqrt(3));
    otherwise
        error('RT:UnknownProto', 'Unknown measurement protocol type `%s`', ptype);
end

end
