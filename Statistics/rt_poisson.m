classdef rt_poisson < rt_statistics
% RT_POISSON The class for working with poisson statistics
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
    methods(Static)
        function k = sample(n, p)
            if numel(p) > 1
                error('RT:PoissDistributionProbCount', 'Only a single result per measurement is allowed for poisson statistics simulation');
            end
            k = poissrnd(p * n);
        end
        
        function f = logL(n, k, p)
            lam = n .* p;
            f = sum(k .* log(lam) - lam);
        end
        
        function df = dlogL(n, k, p)
            df = k ./ p - n;
        end
        
        function mu = logL_mu(~, k)
            mu = sum(k);
        end
        
        function [b, b0] = logL_jmat(n, k, p)
            b = k ./ p - n;
            b0 = sum(n .* p);
        end
        
        function f = chi2(n, k, p)
            ne = p .* n;
            no = k;
            f = sum((ne-no).^2./ne);
        end
        
        function nu = deg_f(clicks)
            nu = length(clicks);
        end
        
        function fi = fisher_information(n, p)
            fi = n ./ p;
        end
    end
end

