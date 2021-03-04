function [dm, c] = rt_pinv(P, freq, r)
%rt_PINV Summary of this function goes here
%   Detailed explanation goes here

dim = size(P, 1);
if nargin < 3
    r = dim;
end

B = rt_meas_matrix(P);
dm = reshape(mldivide(B, freq), dim, dim);
c = rt_purify(dm, r);
dm = c*c';

end

