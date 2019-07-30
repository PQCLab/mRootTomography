function d = rt_chi_theory(chi, M, n)
%rt_CHI_THEORY Summary of this function goes here
%   Detailed explanation goes here

% Chi purification
chi = chi / trace(chi);
r = rank(chi);
[U,D] = svd(chi);
U = U(:,1:r);
D = D(1:r,1:r);
e = U*sqrt(D);

% Information matrix
m = length(M);
if length(n) == 1
    n = ones(1,m)*n;
elseif m ~= length(n)
    error('The length of n should be equal to the number of measurement schemes');
end

s = size(M{1},3);
H = 0;
for j = 1:m
    Hj = 0;
    for i = 1:s
        p = real(trace(M{j}(:,:,i)*chi))+1e-8;
        a = M{j}(:,:,i)*e;
        a = a(:);
        a = [real(a); imag(a)]; 
        Hj = Hj + (a*a')/p;
    end
    H = H + n(j)*Hj;
end
H = 2*H;

% Constraints
R = zeros(size(H,1),s^2);
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

% Finding d
h = svd(U'*H*U);
h(h < 1) = [];
d = 1./(2*s*h);

end

