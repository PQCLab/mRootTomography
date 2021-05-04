function l = rt_lossfun(dm, ex, varargin)
% RT_LOSSFUN Calculates the tomography loss function for a quantum state or a quantum process
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
if isempty(ex.get_field('nshots'))
    ex.set_data('nshots', rt_nshots_divide(1, length(ex.get_field('proto')), 'total'));
end
dF = sum(rt_bound(dm, ex, varargin{:}));
l = dF * sum(ex.get_field('nshots'));
end

