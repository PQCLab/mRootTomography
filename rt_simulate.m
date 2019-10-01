function clicks = rt_simulate(dm, proto, nshots, asymp)
%RT_SIMULATE TODO

if ~iscell(proto)
    proto = num2cell(proto,[1 2]);
end
if nargin < 4
    asymp = false;
end

rt_proto_check(nshots,proto);

m = length(proto);
clicks = cell(1,m);
for j = 1:m
    d = size(proto{j},3);
    p = zeros(d,1);
    for i = 1:d
        p(i) = abs(trace(proto{j}(:,:,i) * dm));
    end
    if asymp
        clicks{j} = nshots(j)*p;
    elseif d == 1
        clicks{j} = poissrnd(nshots(j)*p);
    elseif d == 2
        clicks{j}(1,1) = binornd(nshots(j), p(1));
        clicks{j}(2,1) = nshots(j) - clicks{j}(1,1);
    else
        clicks{j} = mnrnd(nshots(j), p/sum(p))';
    end
end

end

