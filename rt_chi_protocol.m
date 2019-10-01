function [M, P0, M0] = rt_chi_protocol(P0, M0, N)
%rt_CHI_PROTOCOL Measurement protocol for chi-matrix
%   INPUT:
%   P0 - input density matrices
%       string - 'cube' (6 states)
%       OR
%       P0(:,:,i) - density matrix of i-th input state
%   M0 - measurement operators for the state at the output of the gate
%       string - 'pauli' (projective measurement of 3 pauli observables)
%       OR
%       M0{i}(:,:,j) - j-th operator that corresponds to the i-th
%       measurement
%   N - number of qubits (for authomatically generated P0 and M0)
%
%   OUTPUT:
%   M - measurement operators for chi-matrix
%       M{i}(:,:,j) - j-th operator that corresponds to the i-th
%       measurement scheme
%   P0 - input density matrices
%   M0 - measurement operators for the state at the output of the gate

if ischar(P0)
    P0 = rt_proto_preparation(P0, N);
end

if ischar(M0)
    M0 = rt_proto_measurement(M0, N);
end

mp = size(P0,3);
mm = length(M0);

M = cell(1,mp*mm);
for i = 1:mp
    for j = 1:mm        
        M{(i-1)*mm + j} = rt_kron3d(conj(P0(:,:,i)), M0{j});
    end
end

end

