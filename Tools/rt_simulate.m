function clicks = rt_simulate(dm, proto, nshots, stat_type)
%RT_SIMULATE TODO

if nargin < 4
    stat_type = 'auto';
end

d = size(dm, 1);
ex = rt_experiment(d, stat_type, 'state');
ex.set_data('proto', proto, 'nshots', nshots);
clicks = cell(size(ex.proto));
for j = 1:length(ex.proto)
    probs = abs(rt_meas_matrix(ex.proto{j}) * dm(:));
    clicks{j} = ex.sample(probs, ex.nshots(j));
end

end

