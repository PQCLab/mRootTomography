classdef rt_experiment < handle
% RT_EXPERIMENT The class for working with the quantum tomography data
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
    properties
        dim                 % Hilbert space dimension
        obj_type            % Object type (`'state'` or `'process'`)
        stat_type           % Measurements statistics type
        proto = {}          % Measurements operators
        nshots = []         % Measurements repetitions
        clicks = {}         % Number of observed measurements outcomes
        vec_proto = []      % Matrix form of the whole measurements protocol
        vec_nshots = []     % Matrix form of the measurements repetitions
        vec_clicks = []     % Matrix form of the number of observed measurements outcomes
    end
    methods
        function obj = rt_experiment(dim, obj_type, stat_type)
            obj.dim = dim;
            
            if nargin < 3
                stat_type = 'auto';
            end
            switch stat_type
                case {'poly', 'poiss', 'asymp', 'auto'}
                    obj.stat_type = lower(stat_type);
                otherwise
                    error('RT:StatsType', 'Unknown statistics type: `%s`\n Only `poly`, `poiss`, `asymp`, `auto` are available', stat_type);
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
                    % Check whether the process protocol is set as
                    % prepare-measure pairs and transform to chi-matrix
                    % measurement operators
                    if strcmp(obj.obj_type, 'process') && iscell(obj.proto) && all(size(obj.proto) > 1)
                        if size(obj.proto, 2) == 2
                            obj.proto = transpose(obj.proto);
                        end
                        obj.proto = cellfun(@(prep, meas) rt_kron3d(conj(prep), meas), obj.proto(1,:), obj.proto(2,:), 'UniformOutput', false);
                    end
                    % Check whether the protocol is set as 3d array and
                    % transform it ti cell array
                    if ~iscell(obj.proto)
                        dims = size(obj.proto);
                        obj.proto = mat2cell(obj.proto, dims(1), dims(2), ones(1, dims(3)));
                    end
                    obj.proto = reshape(obj.proto, 1, []);
                    % Check whether the statistics type is set to `auto`
                    % and determine statistics according to the protocol
                    % operators
                    if strcmp(obj.stat_type, 'auto')
                        imat = eye(obj.dim);
                        if all(cellfun(@(pr) size(pr, 3) == 1, obj.proto)) % is poiss
                            obj.stat_type = 'poiss';
                        elseif strcmp(obj.obj_type, 'state') && all(cellfun(@(pr) norm(sum(pr, 3) - imat) < 1e-5, obj.proto)) % is state povm
                            obj.stat_type = 'poly';
                        elseif strcmp(obj.obj_type, 'process') && all(cellfun(@(pr) norm(rt_prttrace(sum(pr, 3), [obj.dim, obj.dim], 1) - imat) < 1e-5, obj.proto)) % is process povm
                            obj.stat_type = 'poly';
                        else
                            error('RT:StatsTypeAuto', 'Failed to determine statistics type. Please, specify stat_type manually.');
                        end
                    end
                    continue;
                end
                if strcmp(field, 'clicks')
                    if ~iscell(obj.clicks)
                        obj.clicks = reshape(num2cell(obj.clicks), 1, []);
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
                        obj.nshots = rt_nshots_divide(obj.nshots, size(obj.proto, 2));
                    elseif length(obj.nshots) ~= size(obj.proto, 2)
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
        function p = get_probs_dm(obj, dm, tol)
            if nargin < 3
                tol = 0;
            end
            p = abs(obj.get_field('vec_proto') * dm(:));
            p(p < tol) = tol;
        end
        function p = get_probs_sq(obj, sq, tol)
            if nargin < 3
                tol = 0;
            end
            p = obj.get_probs_dm(sq * sq', tol);
        end
        % ========= Sampling ============
        function clicks = simulate(obj, dm)
            clicks = cell(size(obj.proto));
            for j = 1:numel(obj.proto)
                probs = abs(rt_meas_matrix(obj.proto{j}) * dm(:));
                clicks{j} = obj.sample(probs, obj.nshots(j));
            end
        end
        function k = sample(obj, p, n)
            p = p(:);
            if strcmp(obj.stat_type, 'poly')
                if abs(sum(p) - 1) > 1e-8
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
            elseif strcmp(obj.stat_type, 'poiss')
                k = poissrnd(p * n);
            elseif strcmp(obj.stat_type, 'asymp')
                k = p * n;
            end
            k = k(:);
        end
        % ========= Likelihood ============
        function f = get_logL_dm(obj, dm)
            p = obj.get_probs_dm(dm, 1e-15);
            k = obj.get_field('vec_clicks');
            if strcmp(obj.stat_type, 'poly')
                f = sum(k .* log(p));
            elseif strcmp(obj.stat_type, 'poiss')
                lam = obj.get_field('vec_nshots') .* p;
                f = sum(k .* log(lam) - lam);
            end
        end
        function f = get_logL_sq(obj, sq)
            f = obj.get_logL_dm(sq * sq');
        end
        function df = get_dlogL_sq(obj, sq)
            p = obj.get_probs_sq(sq, 1e-15);
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
                dim2 = obj.dim^2;
                nu = (2*dim2 - obj_rank)*obj_rank - dim2;
            end
            df = df - nu;
        end
    end
end

