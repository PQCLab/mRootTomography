function E = rt_root2kraus(e)

[s2,r] = size(e);
s = sqrt(s2);
E = reshape(e,[s,s,r]);

end

