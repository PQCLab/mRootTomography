classdef rt_poisson_unity < rt_poisson
% RT_POISSON_UNITY The class for working with poisson statistics taking into
% account condition sum_j{n_j * P_j) = const * eye(dim)
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
    methods(Static)        
        function [b, b0] = logL_jmat(~, k, p)
            b = k ./ p;
            b0 = 0;
        end
    end
end

