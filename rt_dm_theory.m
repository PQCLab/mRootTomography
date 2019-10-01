function d = rt_dm_theory(dm, proto, nshots, r)
%RT_DM_THEORY Summary of this function goes here TODO
%   Detailed explanation goes here TODO

if ~iscell(proto)
    proto = num2cell(proto,[1 2]);
end
if nargin < 4 || strcmp(r,'auto')
    r = rank(dm);
end

% Density matrix purification
[U,D] = svd(dm);
U = U(:,1:r);
D = D(1:r,1:r);
D = D / trace(D);
c = U*sqrt(D);

% Information matrix
rt_proto_check(proto, nshots);
[M, n] = rt_data_join(proto, nshots);

H = 0;
pTol = 1e-6;
for j = 1:size(M,3)
    p = real(trace(M(:,:,j)*dm));
    ch = c;
    if p < pTol
        ch = c + (randn(size(c)) + 1j*randn(size(c)))*pTol;
        ch = ch / sqrt(trace(ch'*ch));
        p = real(trace(ch'*M(:,:,j)*ch));
    end
    a = reshape(M(:,:,j)*ch, [], 1);
    a = [real(a); imag(a)]; 
    H = H + n(j)*(a*a')/p;
end
H = 2*H;

% Finding d
h = svd(H);
h(1) = [];
ind0 = find(h < 1e-4);
if length(ind0) > r^2
   warning('Information matrix has more than r^2 zero eigenvalues. Extra zeros are ignored');
end
h(ind0) = [];
d = 1./(2*h);

end

