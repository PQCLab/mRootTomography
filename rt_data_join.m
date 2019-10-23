function [M, n, k] = rt_data_join(proto, nshots, clicks)

M = [];
n = [];
k = [];
for j = 1:length(proto)
    M = cat(3, M, proto{j});
    if nargin > 1
        n = vertcat(n, repmat(nshots(j),size(proto{j},3),1));
    end
    if nargin > 2
        k = vertcat(k, clicks{j});
    end
end

end

