classdef rt_experiment
    %RT_STATS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dim
        stat_type
        proto = {}
        nshots = []
        clicks = {}
        vec_proto = []
        vec_nshots = []
        vec_clicks = []
    end
    
    methods
        function obj = rt_experiment(dim, stat_type)
            obj.dim = dim;
            switch stat_type
                case {'poly', 'poiss'}
                    obj.stat_type = stat_type;
                otherwise
                    error('RT:StatsType', 'Unknown statistics type: `%s`\n Only `poly` and `poiss` are available', stat_type);
            end
        end
        
        function obj = set_data(obj, varargin)
            for j = 1:2:length(varargin)
                obj.(varargin{j}) = varargin{j+1};
            end
        end
        
        function data = get_field(obj, field)
            if strcmp(field, 'vec_proto') && isempty(obj.vec_proto)
                obj.vec_proto = rt_meas_matrix(cat(3, obj.proto{:}));
            elseif strcmp(field, 'vec_nshots') && isempty(obj.vec_nshots)
                n = arrayfun(@(nj,Mj) ones(size(Mj{1},3),1)*nj, obj.nshots(:), obj.proto(:), 'UniformOutput', false);
                obj.vec_nshots = vertcat(n{:});
            elseif strcmp(field, 'vec_clicks') && isempty(obj.vec_clicks)
                k = cellfun(@(kj) kj(:), obj.clicks, 'UniformOutput', false);
                obj.vec_clicks = vertcat(k{:});
            end
            data = obj.(field);
        end
        
        function p = get_probs_dm(obj, dm)
            p = abs(obj.get_field('vec_proto') * dm(:));
            p(p < 1e-15) = 1e-15;
        end
        
        function p = get_probs_sq(obj, sq)
            p = obj.get_probs_dm(sq * sq');
        end
        
        
        % ========= Sampling ============
        function k = sample(obj, p, n)
            p = p(:);
            if strcmpi(obj.stat_type, 'poly')
                if n > 1e5 % normal approximation for performance
                    mu = p*n;
                    sigma = (-p*p' + diag(p))*n;
                    k = vec(round(mvnrnd(mu, sigma)));
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
                        k = vec(mnrnd(n, p));
                    end
                end
            elseif strcmpi(obj.stat_type, 'poly')
                k = poissrnd(p*n);
            end
        end
        
        
        % ========= Likelihood ============
        function f = get_logL_dm(obj, dm)
            p = obj.get_probs_dm(dm);
            k = obj.get_field('vec_clicks');
            if strcmpi(obj.stat_type, 'poly')
                f = sum(k.*log(p));
            elseif strcmpi(obj.stat_type, 'poiss')
                lam = obj.get_field('vec_nshots') .* p;
                f = sum(k.*log(lam) - lam);
            end
        end
        
        function f = get_logL_sq(obj, sq)
            f = obj.get_logL_dm(sq * sq');
        end
        
        function df = get_dlogL_sq(obj, sq)
            p = obj.get_probs_sq(sq);
            k = obj.get_field('vec_clicks');
            B = obj.get_field('vec_proto');
            a = k ./ p;
            if strcmpi(obj.stat_type, 'poly')
                J = reshape(B' * a, size(sq, 1), []);
                df = 2 * J * sq;
            elseif strcmpi(obj.stat_type, 'poiss')
                n = obj.get_field('vec_nshots');
                J = reshape(B' * (a - n), size(sq, 1), []);
                df = 2 * J * sq;
            end
        end
        
        % ========= Chi-squared ============
    end
end

