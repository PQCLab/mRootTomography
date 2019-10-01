function M = rt_proto_measurement(type, N, separate)
%rt_PROTOCOL_MEASUREMENT Generates measurement protocol for the state
%   INPUT:
%   type - protocol type
%       'pauli' - Pauli measurements (measurement of s_x, s_y, s_z)
%       'tetra' - Measurement of tetrahedron observables
%   N - number of qubits (default: 1)
%   separate - true if each measurement operator is measured separately (default: false)
%
%   OUTPUT:
%   M - measurement operators
%       M0{i}(:,:,j) - j-th operator that corresponds to the i-th
%       measurement

if nargin < 2
    N = 1;
end
if nargin < 3
    separate = false;
end

switch type
    case 'pauli'
        M = cell(1,3);
        M{1}(:,:,1) = [1 1; 1 1]/2;
        M{1}(:,:,2) = [1 -1; -1 1]/2;
        M{2}(:,:,1) = [1 -1j; 1j 1]/2;
        M{2}(:,:,2) = [1 1j; -1j 1]/2;
        M{3}(:,:,1) = [1 0; 0 0];
        M{3}(:,:,2) = [0 0; 0 1];
    case {'tetrahedron', 'octahedron'}
        M = cell(1,4);
        M{1}(:,:,1) = [sqrt(3)+1 sqrt(2)*exp(1j*1*pi/4); sqrt(2)*exp(-1j*1*pi/4) sqrt(3)-1]/(2*sqrt(3));
        M{1}(:,:,2) = [sqrt(3)-1 sqrt(2)*exp(1j*5*pi/4); sqrt(2)*exp(-1j*5*pi/4) sqrt(3)+1]/(2*sqrt(3));
        M{2}(:,:,1) = [sqrt(3)-1 sqrt(2)*exp(1j*7*pi/4); sqrt(2)*exp(-1j*7*pi/4) sqrt(3)+1]/(2*sqrt(3));
        M{2}(:,:,2) = [sqrt(3)+1 sqrt(2)*exp(1j*3*pi/4); sqrt(2)*exp(-1j*3*pi/4) sqrt(3)-1]/(2*sqrt(3));
        M{3}(:,:,1) = [sqrt(3)+1 sqrt(2)*exp(1j*5*pi/4); sqrt(2)*exp(-1j*5*pi/4) sqrt(3)-1]/(2*sqrt(3));
        M{3}(:,:,2) = [sqrt(3)-1 sqrt(2)*exp(1j*1*pi/4); sqrt(2)*exp(-1j*1*pi/4) sqrt(3)+1]/(2*sqrt(3));
        M{4}(:,:,1) = [sqrt(3)-1 sqrt(2)*exp(1j*3*pi/4); sqrt(2)*exp(-1j*3*pi/4) sqrt(3)+1]/(2*sqrt(3));
        M{4}(:,:,2) = [sqrt(3)+1 sqrt(2)*exp(1j*7*pi/4); sqrt(2)*exp(-1j*7*pi/4) sqrt(3)-1]/(2*sqrt(3));
    otherwise
        error('Unknown measurement protocol');
end

if separate
    M = cat(3, M{:});
    if strcmp(type, 'tetrahedron')
        M = M(:,:,1:2:end);
    end
    M = num2cell(M,[1 2]);
end

M1 = M;
m1 = length(M1);
for i = 2:N
    m = length(M);
    Mnew = cell(1,m*m1);
    for j = 1:m
        for k = 1:m1
            Mnew{(j-1)*m1 + k} = rt_kron3d(M{j}, M1{k});
        end
    end
    M = Mnew;
end

end
