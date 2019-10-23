function P = rt_proto_preparation(type, N)
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

P = 1;
for j = 1:N
    P = rt_kron3d(P, P0);
end

end

