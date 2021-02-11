classdef rt_experiment < handle
    %RT_STATS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dim
        stat_type
        obj_type
        proto = {}
        nshots = []
        clicks = {}
        vec_proto = []
        vec_nshots = []
        vec_clicks = []
    end
    
    methods
        function obj = rt_experiment(dim, stat_type, obj_type)
            obj.dim = dim;
            
            if nargin < 2
                stat_type = 'poiss';
            end
            switch stat_type
                case {'poly', 'poiss', 'bino', 'asymp', 'auto'}
                    obj.stat_type = lower(stat_type);
                otherwise
                    error('RT:StatsType', 'Unknown statistics type: `%s`\n Only `poly`, `poiss`, `bino`, `asymp`, `auto` are available', stat_type);
            end
            
            if nargin < 3
                obj_type = 'state';
            end
            switch obj_type
                case {'state', 'process'}
                    obj.obj_type = lower(obj_type);
                otherwise
                    error('RT:ObjectType', 'Unknown object type: `%s`\n Only `state`, `process` are available', obj_type);
            end
        end
        
        function obj = set_data(obj, varargin)
            for j = 1:2:length(varargin)
                field = lower(varargin{j});
                obj.(field) = varargin{j+1};
                if strcmp(field, 'proto')
                    if ~iscell(obj.proto)
                        obj.proto = reshape(mat2cell(obj.proto, obj.dim, obj.dim, ones(1, size(obj.proto, 3))), 1, []);
                    end
                    if strcmp(obj.stat_type, 'auto')
                        imat = eye(obj.dim);
                        if all(cellfun(@(pr) norm(sum(pr, 3) - imat) < 1e-8, obj.proto)) % is povm
                            obj.stat_type = 'poly';
                        elseif all(cellfun(@(pr) size(pr, 3) == 1, obj.proto)) % is poiss
                            obj.stat_type = 'poiss';
                        else
                            error('RT:StatsTypeAuto', 'Failed to determine statistics type');
                        end
                    end
                    continue;
                end
                if strcmp(field, 'clicks')
                    if ~iscell(obj.clicks)
                        obj.clicks = num2cell(obj.clicks);
                    end
                    continue;
                end
            end
            if ~isempty(obj.nshots)
                if strcmp(obj.stat_type, 'auto') && any(isinf(obj.nshots))
                    obj.nshots = ones(size(obj.nshots));
                    obj.stat_type = 'asymp';
                end
                if ~isempty(obj.proto)
                    if length(obj.nshots) == 1
                        obj.nshots = rt_nshots_devide(obj.nshots, length(obj.proto));
                    elseif length(obj.nshots) ~= length(obj.proto)
                        error('RT:ExpNumberMismatch', 'Length of nshots array does not match length of proto array');
                    end
                end
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
            if strcmp(obj.stat_type, 'poly')
                p = p / sum(p);
                if n > 1e5 % normal approximation for performance
                    mu = p*n;
                    sigma = (-p*p' + diag(p))*n;
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
            elseif strcmp(obj.stat_type, 'poiss')
                k = poissrnd(p * n);
            elseif strcmp(obj.stat_type, 'bino')
                k = binornd(n, p);
            elseif strcmp(obj.stat_type, 'asymp')
                k = p * n;
            end
            k = k(:);
        end
        
        
        % ========= Likelihood ============
        function f = get_logL_dm(obj, dm)
            p = obj.get_probs_dm(dm);
            k = obj.get_field('vec_clicks');
            if strcmp(obj.stat_type, 'poly')
                f = sum(k .* log(p));
            elseif strcmp(obj.stat_type, 'poiss')
                lam = obj.get_field('vec_nshots') .* p;
                f = sum(k .* log(lam) - lam);
            elseif strcmp(obj.stat_type, 'bino')
                n = obj.get_field('vec_nshots');
                f = sum(k .* log(p) + (n - k) .* log(1 - p));
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
            elseif strcmpi(obj.stat_type, 'bino')
                n = obj.get_field('vec_nshots');
                J = reshape(B' * (a - (n - k)./(1 - p)), size(sq, 1), []);
                df = 2 * J * sq;
            end
        end
        
        % ========= Chi-squared ============
        function f = get_chi2_dm(obj, dm)
            n_expected = obj.get_probs_dm(dm) .* obj.get_field('vec_nshots');
            n_observed = obj.get_field('vec_clicks');
            f = sum((n_expected-n_observed).^2./n_expected);
        end
        
        function df = get_df(obj, obj_rank)
            df = length(obj.get_field('vec_clicks'));
            if strcmpi(obj.stat_type, 'poly')
                df = df - length(obj.get_field('clicks'));
            end
            if strcmpi(obj.obj_type, 'state')
                nu = (2*obj.dim - obj_rank)*obj_rank - 1;
            elseif strcmpi(obj.obj_type, 'process')
                nu = (2*obj.dim - obj_rank)*obj_rank - obj.dim;
            end
            df = df - nu;
        end
    end
end

