classdef rt_binomial < rt_statistics
% RT_BINOMIAL The class for working with binomial statistics
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
    methods(Static)
        function k = sample(n, p)
            if numel(p) > 1
                error('RT:BinoDistributionProbCount', 'Only a single result per measurement is allowed for binomial statistics simulation');
            end
            k = binornd(n, p);
        end
        
        function f = logL(n, k, p)
            f = sum(k .* log(p) + (n - k) .* log(1 - p));
        end
        
        function df = dlogL(n, k, p)
            df = k ./ p - (n - k) ./ (1 - p);
        end
        
        function mu = logL_mu(~, k)
            mu = sum(k);
        end
        
        function [b, b0] = logL_jmat(n, k, p)
            b = k ./ p - (n - k) ./ (1 - p);
            b0 = sum((n - k) .* p ./ (1 - p));
        end
        
        function f = chi2(n, k, p)
            ne = [p .* n; (1 - p) .* n];
            no = [k; n - k];
            f = sum((ne-no).^2./ne);
        end
        
        function nu = deg_f(clicks)
            nu = length(clicks);
        end
        
        function fi = fisher_information(n, p)
            fi = n ./ p ./ (1 - p);
        end
    end
end

