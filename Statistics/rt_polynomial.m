classdef rt_polynomial < rt_statistics
% RT_POLYNOMIAL The class for working with polynomial statistics
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
    methods(Static)
        function k = sample(n, p)
            p = p(:);
            if abs(sum(p) - 1) > 1e-5
                error('RT:PolyDistributionNorm', 'For simulating polynomial statistics probabilities in each measurement should sum to unity');
            end
            p = p / sum(p);
            if n > 1e5 % normal approximation for performance
                mu = p*n;
                sigma = (-p*p' + diag(p))*n;
                [u, d] = eigs(sigma, rank(sigma));
                d(d < 0) = 0;
                sigma = u * d * u';
                k = reshape(round(mvnrnd(mu, sigma)), [], 1);
                k(k < 0) = 0;
                if sum(k) > n
                    kmind = find(k == max(k),1);
                    k(kmind) = k(kmind) - (sum(k) - n);
                else
                    k(end) = n - sum(k(1:(end-1)));
                end
            else
                if length(p) == 2
                    k = zeros(2, 1);
                    k(1) = binornd(n, p(1));
                    k(2) = n - k(1);
                else
                    k = mnrnd(n, p);
                end
            end
        end
        
        function f = logL(~, k, p)
            f = sum(k .* log(p));
        end
        
        function df = dlogL(~, k, p)
            df = k ./ p;
        end
        
        function mu = logL_mu(~, k)
            mu = sum(k);
        end
        
        function [b, b0] = logL_jmat(~, k, p)
            b = k ./ p;
            b0 = 0;
        end
        
        function f = chi2(n, k, p)
            ne = p .* n;
            no = k;
            f = sum((ne-no).^2./ne);
        end
        
        function nu = deg_f(clicks)
            nu = sum(cellfun(@(k) length(k)-1, clicks));
        end
        
        function fi = fisher_information(n, p)
            fi = n ./ p;
        end
    end
end

