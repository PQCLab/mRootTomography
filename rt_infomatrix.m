function H = rt_infomatrix(dm, proto, nshots, varargin)
% RT_INFOMATRIX TODO
p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'dm');
addRequired(p, 'proto');
addRequired(p, 'nshots');
addOptional(p, 'rank', 'dm');
parse(p,dm,proto,nshots,varargin{:});
opt = p.Results;

if ischar(opt.rank) && strcmp(opt.rank,'dm')
    opt.rank = rank(dm);
end
proto = rt_proto_check(proto, nshots);

c = rt_purify(dm,opt.rank);
[M, n] = rt_data_join(proto, nshots);
B = rt_meas_matrix(M);

pTol = 1e-10;
Ntries = 100;
for i = 1:Ntries
    prob = abs(B*reshape(c*c',[],1));
    if any(prob < pTol)
        if i == Ntries
            warning('Failed to find non-singular state');
        else
            c = c + (randn(size(c)) + 1j*randn(size(c)))*sqrt(pTol);
            c = c / sqrt(trace(c'*c));
        end
    else
        break;
    end
end

H = 0;
for j = 1:size(M,3)
    a = reshape(M(:,:,j)*c, [], 1);
    a = [real(a); imag(a)]; 
    H = H + n(j)*(a*a')/prob(j);
end
H = 2*H;

end

