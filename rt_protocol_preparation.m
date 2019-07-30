function P = rt_protocol_preparation(type, N, Ep)
%rt_PROTOCOL_PREPARATION Generates preparation protocol
%   INPUT:
%   type - protocol type
%       'cube' - Cube protocol (pure eigenstates of s_x, s_y, s_z)
%   N - number of qubits
%   Ep (optional) - NG-model fuzzy measurement operators (see paper)
%       Ep(:,:,k) - k-th Kraus operator
%
%   OUTPUT:
%   P - preparation density matrices
%       P0(:,:,i) - density matrix of the i-th input state

if strcmpi(type, 'cube')
    P0 = zeros(2,2,6);
    P0(:,:,1) = [1 1; 1 1]/2;
    P0(:,:,2) = [1 -1; -1 1]/2;
    P0(:,:,3) = [1 -1j; 1j 1]/2;
    P0(:,:,4) = [1 1j; -1j 1]/2;
    P0(:,:,5) = [1 0; 0 0];
    P0(:,:,6) = [0 0; 0 1];
else
    error('Unknown preparation protocol');
end

if nargin < 3
    Ep = cell(1,N);
    [Ep{:}] = deal(eye(2));
end

P = 1;
for j = 1:N
    Ptmp = P0;
    for k = 1:size(Ptmp,3)
        Ptmp(:,:,k) = rt_kraus(Ptmp(:,:,k), Ep{j});
    end
    P = superkron(P, Ptmp);
end

end

