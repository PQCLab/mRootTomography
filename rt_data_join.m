function [M, n, k] = rt_data_join(proto, nshots, clicks)

M = [];
n = [];
k = [];
for j = 1:length(proto)
    M = cat(3, M, proto{j});
    n = vertcat(n, repmat(nshots(j),size(proto{j},3),1));
    if nargin == 3
        k = vertcat(k, clicks{j});
    end
end

end

