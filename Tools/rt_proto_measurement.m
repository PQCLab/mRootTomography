function proto = rt_proto_measurement(ptype, varargin)
% RT_PROTO_MEASUREMENT Generates the measurements operators of a specific type
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
op.dim = nan;
op.modifier = 'none';
op.nsub = 1;
for ja = 1:2:length(varargin)
    op.(lower(varargin{ja})) = varargin{ja + 1};
end

switch ptype
    case 'mub'
        pname = ['mub', num2str(op.dim)];
        mubs = load('mubs.mat', pname);
        mubs = mubs.(pname);
        proto = cell(1, op.dim + 1);
        for j = 1:(op.dim+1)
            proto{j} = zeros(op.dim, op.dim, op.dim);
            for k = 1:op.dim
                proto{j}(:,:,k) = mubs(:,k,j) * mubs(:,k,j)';
            end
        end
        if strcmpi(op.modifier, 'operator')
            proto = cat(3, proto{:});
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
        if strcmpi(op.modifier, 'operator')
            proto = cat(3, proto{:});
        elseif strcmpi(op.modifier, 'operator+')
            proto = cellfun(@(pr) pr(:,:,1), proto, 'UniformOutput', false);
            proto = cat(3, proto{:});
        elseif strcmpi(op.modifier, 'operator-')
            proto = cellfun(@(pr) pr(:,:,2), proto, 'UniformOutput', false);
            proto = cat(3, proto{:});
        end
    otherwise
        error('RT:UnknownProto', 'Unknown measurement protocol type `%s`', ptype);
end

if ~iscell(proto)
    dims = size(proto);
    proto = reshape(mat2cell(proto, dims(1), dims(2), ones(1, dims(3))), 1, []);
end

if op.nsub > 1
    proto0 = transpose(proto);
    for jn = 2:op.nsub
        proto = repmat(proto, length(proto0), 1);
        proto0r = repmat(proto0, 1, length(proto));
        for rp = 1:numel(proto)
            proto{rp} = rt_kron3d(proto{rp}, proto0r{rp});
        end
        proto = reshape(proto, 1, []);
    end
end

end
