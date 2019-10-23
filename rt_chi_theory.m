function d = rt_chi_theory(chi, proto, nshots, r)
%rt_CHI_THEORY TODO

chi = chi / trace(chi);
if nargin < 4 || strcmp(r,'auto')
    r = rank(chi);
end
proto = rt_proto_check(proto, nshots);
H = rt_infomatrix(chi, proto, nshots, r);

% Constraints
e = rt_purify(dm,r);
s2 = size(chi,1);
s = sqrt(s2);
R = zeros(size(H,1),s2);
ind = @(i,j) (i-1)*s + j;
r2 = 1:s;
for m1 = 1:s
    for n1 = 1:s
        for r1 = 1:s
            indr = ind(r1,r2);
            [deriv_re, deriv_im] = deal(zeros(size(e)));
            if r1 == n1
                deriv_re(indr,:) = deriv_re(indr,:) + e(ind(m1,r2),:);
                deriv_im(indr,:) = deriv_im(indr,:) - e(ind(m1,r2),:)*1j;
            end    
            if r1 == m1
                deriv_re(indr,:) = deriv_re(indr,:) + conj(e(ind(n1,r2),:));
                deriv_im(indr,:) = deriv_im(indr,:) + conj(e(ind(n1,r2),:))*1j;
            end
            R(:,ind(m1,n1)) = R(:,ind(m1,n1)) + [deriv_re(:); deriv_im(:)];
        end
    end
end
U = null(R');
H = U'*H*U;

% Finding d
h = svd(H);
ind0 = find(h < 1e-4);
if length(ind0) > r^2
   warning('Information matrix has more than r^2 zero eigenvalues. Extra zeros are ignored');
end
h(ind0) = [];
d = 1./(2*s*h);

end

