function M = rt_protocol_measurement(type, N, Em)
%rt_PROTOCOL_MEASUREMENT Generates measurement protocol for the state
%   INPUT:
%   type - protocol type
%       'pauli' - Pauli measurements (measurement of s_x, s_y, s_z)
%       'tetra' - Measurement of tetrahedron observables
%   N - number of qubits
%   Em (optional) - GN-model fuzzy measurement operators (see paper)
%       Em(:,:,k) - k-th Kraus operator
%
%   OUTPUT:
%   M - measurement operators
%       M0{i}(:,:,j) - j-th operator that corresponds to the i-th
%       measurement

if nargin < 2
    N = 1;
end

if ischar(type)
    switch type
        case 'pauli'
            Mm = cell(1,3);
            Mm{1}(:,:,1) = [1 1; 1 1]/2;
            Mm{1}(:,:,2) = [1 -1; -1 1]/2;
            Mm{2}(:,:,1) = [1 -1j; 1j 1]/2;
            Mm{2}(:,:,2) = [1 1j; -1j 1]/2;
            Mm{3}(:,:,1) = [1 0; 0 0];
            Mm{3}(:,:,2) = [0 0; 0 1];
        case 'tetra'
            Mm = cell(1,4);
            Mm{1}(:,:,1) = [sqrt(3)+1 sqrt(2)*exp(1j*1*pi/4); sqrt(2)*exp(-1j*1*pi/4) sqrt(3)-1]/(2*sqrt(3));
            Mm{1}(:,:,2) = [sqrt(3)-1 sqrt(2)*exp(1j*5*pi/4); sqrt(2)*exp(-1j*5*pi/4) sqrt(3)+1]/(2*sqrt(3));
            Mm{2}(:,:,1) = [sqrt(3)-1 sqrt(2)*exp(1j*7*pi/4); sqrt(2)*exp(-1j*7*pi/4) sqrt(3)+1]/(2*sqrt(3));
            Mm{2}(:,:,2) = [sqrt(3)+1 sqrt(2)*exp(1j*3*pi/4); sqrt(2)*exp(-1j*3*pi/4) sqrt(3)-1]/(2*sqrt(3));
            Mm{3}(:,:,1) = [sqrt(3)+1 sqrt(2)*exp(1j*5*pi/4); sqrt(2)*exp(-1j*5*pi/4) sqrt(3)-1]/(2*sqrt(3));
            Mm{3}(:,:,2) = [sqrt(3)-1 sqrt(2)*exp(1j*1*pi/4); sqrt(2)*exp(-1j*1*pi/4) sqrt(3)+1]/(2*sqrt(3));
            Mm{4}(:,:,1) = [sqrt(3)-1 sqrt(2)*exp(1j*3*pi/4); sqrt(2)*exp(-1j*3*pi/4) sqrt(3)+1]/(2*sqrt(3));
            Mm{4}(:,:,2) = [sqrt(3)+1 sqrt(2)*exp(1j*7*pi/4); sqrt(2)*exp(-1j*7*pi/4) sqrt(3)-1]/(2*sqrt(3));
        otherwise
            error('Unknown measurement protocol');
    end
else
    Mm = type;
end

if nargin == 3
    if N > 1 && length(Em) == 1
        Etmp = Em;
        Em = cell(1,N);
        [Em{:}] = deal(Etmp);
    end
else
    Em = cell(1,N);
    [Em{:}] = deal(eye(2));
end

mm = length(Mm);
M = noisify_M(Mm, Em{1});
for i = 2:N
    m = length(M);
    Mn = cell(1,m*mm);
    Mmn = noisify_M(Mm,Em{i});
    for j = 1:m
        for k = 1:mm
            Mn{(j-1)*mm + k} = superkron(M{j}, Mmn{k});
        end
    end
    M = Mn;
end

end

function M = noisify_M(M, E)

mm = length(M);
re = size(E, 3);
for j = 1:mm
    for i = 1:2
        Mtmp = 0;
        for k = 1:re
            Mtmp = Mtmp + E(:,:,k)'*M{j}(:,:,i)*E(:,:,k);
        end
        M{j}(:,:,i) = Mtmp;
    end
end

end

