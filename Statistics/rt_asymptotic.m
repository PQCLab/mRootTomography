classdef rt_asymptotic < rt_poisson
% RT_ASYMPTOTIC The class for working with asymptotic statistics
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
    methods(Static)
        function k = sample(n, p)
            k = p * n;
        end
    end
end