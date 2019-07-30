function U = rt_randunitary(s)
%RANDUNITARY Generates random unitary matrix of dimension sxs

if nargin < 1
    s = 2;
end

[U,~] = svd(rand(s) + 1j*rand(s));

end

