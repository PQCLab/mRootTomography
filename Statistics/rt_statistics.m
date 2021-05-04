classdef rt_statistics < handle
% RT_STATISTICS The abstract class for working with different statistics types
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
    properties (Constant)
        buildin = struct( ...
            'polynomial',       @rt_polynomial, ...
            'binomial',         @rt_binomial, ...
            'poisson',          @rt_poisson, ...
            'poisson_unity',    @rt_poisson_unity, ...
            'asymptotic',       @rt_asymptotic ...
        )
    end
    
    methods(Static)
        function obj = get(name)
            obj = rt_statistics.buildin.(name)();
        end
    end
    
    methods(Abstract)
        k = sample(n, p)
        f = logL(n, k, p)
        df = dlogL(n, k, p)
        mu = logL_mu(n, k)
        [b, b0] = logL_jmat(n, k, p)
        f = chi2(n, k, p)
        nu = deg_f(clicks)
        fi = fisher_information(n, p)
    end
end

