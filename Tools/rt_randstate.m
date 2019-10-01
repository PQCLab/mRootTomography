function state = randstate(s, mode, r)
%RANDSTATE Generates random state
%   Generates random state of specific dimention s

if nargin < 2
    mode = 'complex';
end
if nargin < 1
    s = 2;
end

p = [];
if length(s) > 1
    p = s;
    s = length(p);
    mode = 'mixed';
end

if strcmp(mode, 'mixed')
    if nargin < 3
        r = s;
    end
    
    [U,~,~] = svd(rand(s)+1j*rand(s));
    
    if isempty(p)
        p = rand(1,r);
        p = p/sum(p);
        p = [p zeros(1,s-length(p))];
    end
    
    state = U*diag(p)*U';
else
    state = rand(s,1);
    if strcmp(mode, 'complex')
        state = state + 1j*rand(s,1);
    end
    state = state / sqrt(state'*state);
end

end
