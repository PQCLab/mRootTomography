function clicks = rt_simulate(dm, proto, nshots, asymp)
%RT_SIMULATE TODO

[proto,nshots] = rt_proto_check(proto,nshots);

if nargin < 4
    asymp = false;
end
if nargin < 3
    nshots = rt_nshots_devide(1,length(proto),'equal');
    asymp = true;
end


m = length(proto);
clicks = cell(1,m);
for j = 1:m
    d = size(proto{j},3);
    p = zeros(d,1);
    for i = 1:d
        p(i) = real(trace(proto{j}(:,:,i) * dm));
    end
    if asymp
        clicks{j} = nshots(j)*p;
    elseif d == 1
        clicks{j} = poissrnd(nshots(j)*p);
    elseif d == 2
        clicks{j}(1,1) = binornd(nshots(j), p(1)/sum(p));
        clicks{j}(2,1) = nshots(j) - clicks{j}(1,1);
    else
        clicks{j} = mnrnd(nshots(j), p/sum(p))';
    end
end

end

