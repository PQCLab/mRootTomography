function [dm, c] = rt_pinv(M, p, r)
%rt_PINV Summary of this function goes here
%   Detailed explanation goes here

s = size(M,1);
B = rt_meas_matrix(M);

dm = reshape(B\p, s, s);
c = rt_purify(dm, r);
dm = c*c';

end

