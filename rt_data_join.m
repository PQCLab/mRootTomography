function [M, n, k] = rt_data_join(proto, nshots, clicks)

M = cat(3, proto{:});
if nargin > 1
    n = arrayfun(@(nj,Mj) ones(size(Mj{1},3),1)*nj, nshots(:), proto(:), 'UniformOutput', false);
    n = vertcat(n{:});
end
if nargin > 2
    k = cellfun(@(kj) kj(:), clicks, 'UniformOutput', false);
    k = vertcat(k{:});
end

end

