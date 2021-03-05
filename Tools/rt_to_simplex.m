function p = rt_to_simplex(p)
% RT_TO_SIMPLEX Calculates the projection of a real-valued vector onto standard simplex
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
d = length(p);
[p, srt_idx] = sort(p, 'descend');
pcum = cumsum(p);
mu = (pcum - 1) ./ vec(1:d);
idx = find(p - mu > 0, 1, 'last');
p = max(p - mu(idx), 0);
p(srt_idx) = p;
end

