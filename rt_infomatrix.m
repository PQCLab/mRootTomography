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

end

