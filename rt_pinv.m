function [dm, c] = rt_pinv(M, p, r)
%rt_PINV Summary of this function goes here
%   Detailed explanation goes here

s = size(M,1);
B = rt_meas_matrix(M);

% warning('off','MATLAB:rankDeficientMatrix');
% warning('off','MATLAB:nearlySingularMatrix');
dm = reshape(pinv(B)*p, s, s);
% warning('on','MATLAB:rankDeficientMatrix');
% warning('on','MATLAB:nearlySingularMatrix');
c = rt_purify(dm, r);
dm = c*c';

end

