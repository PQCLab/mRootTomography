function d = rt_chi_theory(chi, proto, nshots, r)
%rt_CHI_THEORY TODO

if ~iscell(proto)
    proto = num2cell(proto,[1 2]);
end
if nargin < 4 || strcmp(r,'auto')
    r = rank(chi);
end

% Chi purification
chi = chi / trace(chi);
[U,D] = svd(chi);
U = U(:,1:r);
D = D(1:r,1:r);
c = U*sqrt(D);

% Information matrix
rt_proto_check(proto, nshots);
[M, n] = rt_data_join(proto, nshots);

H = 0;
pTol = 1e-6;
for j = 1:size(M,3)
    p = real(trace(M(:,:,j)*chi));
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

% Constraints
s2 = size(chi,1);
s = sqrt(s2);
R = zeros(size(H,1),s2);
ind = @(i,j) (i-1)*s + j;
r2 = 1:s;
for m1 = 1:s
    for n1 = 1:s
        for r1 = 1:s
            indr = ind(r1,r2);
            [deriv_re, deriv_im] = deal(zeros(size(c)));
            if r1 == n1
                deriv_re(indr,:) = deriv_re(indr,:) + c(ind(m1,r2),:);
                deriv_im(indr,:) = deriv_im(indr,:) - c(ind(m1,r2),:)*1j;
            end    
            if r1 == m1
                deriv_re(indr,:) = deriv_re(indr,:) + conj(c(ind(n1,r2),:));
                deriv_im(indr,:) = deriv_im(indr,:) + conj(c(ind(n1,r2),:))*1j;
            end
            R(:,ind(m1,n1)) = R(:,ind(m1,n1)) + [deriv_re(:); deriv_im(:)];
        end
    end
end
U = null(R');

% Finding d
h = svd(U'*H*U);
h(h < 1) = [];
d = 1./(2*s*h);

end

