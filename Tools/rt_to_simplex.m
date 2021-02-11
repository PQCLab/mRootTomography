function p = rt_to_simplex(p)
%rt_to_simplex Project real non-positive probability distribution to
%positive one

d = length(p);
b = sum(p);
[p, srt_idx] = sort(p, 'descend');
pcum = cumsum(p);
mu = (pcum - b) ./ vec(1:d);
idx = find(p - mu > 0, 1, 'last');
p = max(p - mu(idx), 0);
p(srt_idx) = p;

end

