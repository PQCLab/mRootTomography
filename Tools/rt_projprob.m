function p = rt_projprob(p)
%PROJPROB Project real non-positive probability distribution no positive

a = 0;
[ps,ind] = sort(p,'descend');
for i = length(ps):-1:1
    if ps(i) + a/i >= 0
        ps(1:i) = ps(1:i)+a/i;
        break;
    end
    a = a + ps(i);
    ps(i) = 0;
end
p(ind) = ps;

end

