function p = rt_projprob(p)
%PROJPROB Project real non-positive probability distribution no positive

for i = 1:1000
    in = find(p<0);
    mn = length(in);
    if mn == 0
        break;
    end
    
    ip = find(p>0);
    mp = length(ip);
    
    dp = -sum(p(in))/mp;
    p(in) = 0;
    p(ip) = p(ip)-dp;
end

if ~isempty(find(p<0, 1))
    error('Failed to project probability distribution');
end

end

